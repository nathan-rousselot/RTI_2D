import numpy as np
import matplotlib.pyplot as plt
import yt
import os
import sys
from scipy.constants import m_p, k, R

# plt.rcParams['text.usetex'] = True # Use latex formatting in matplotlib figures 

def getParfileValues(Atwoods: float, Reynolds: float, lamba: float = 0.4, rholight: float = 1.0, taumax: float = 6, he_abundance = 0.1):
    # Compute parameters for in the par file
    # Inputs:
    #   Atwoods number
    #   Reynolds number
    #   lambda = perturbation wavelength
    #   rholight = density of light fluid
    #   taumax = simulation length in dimensionless units
    # Outputs:
    #   vc_mu: viscosity related in &vc_list
    #   time_max: simulation length in physical units in &stoplist

    rhodens=rholight*(1.0+Atwoods)/(1.0-Atwoods)
    vc_mu = lamba*np.sqrt(Atwoods*lamba/(Atwoods+1.0))/(2*Reynolds/(rholight + rhodens))
    time_max = taumax/np.sqrt(Atwoods/lamba)

    unit_time=np.sqrt(lamba/Atwoods)
    unit_velocity=np.sqrt(lamba*Atwoods/(1.0+Atwoods))
    unit_numberdensity=1e24
    unit_density = (1 + 4*he_abundance)*m_p*unit_numberdensity
    unit_length = unit_velocity*unit_time
    unit_mass = unit_density*np.power(unit_length, 3)
    unit_pressure = unit_density*np.power(unit_velocity, 2)
    unit_temperature = unit_pressure/((2 + 3*he_abundance)*unit_numberdensity*k)
    unit_energy = unit_mass*np.power(unit_length, 2)/np.power(unit_time, 2)

    Mp = 1.00727646627*1e-3 # proton molar mass

    vc_mu = vc_mu*unit_length*unit_time/unit_mass

    tc_k_para = vc_mu*5/2*R/Mp # with units
    tc_k_para = tc_k_para*unit_temperature*unit_mass/unit_energy # make unitless

    print(f'{unit_time=}')
    print(f'{unit_velocity=}')
    print(f'{unit_numberdensity=}')
    print(f'{vc_mu=}')
    print(f'{tc_k_para=}')
    print(f'{time_max=}')

    return vc_mu, tc_k_para, time_max

getParfileValues(0.04, 1000)

def getBSHeigthVelocity(fileName: str, level: int = 2):
    # Get the Bubble and spike height and velocity
    # Inputs:
    #   fileName: str: name of .dat file
    #   level: int: level at which to interpolate
    # Outputs:
    #   Buuble/spike heights/velocities: double: hB, FrB, hS, FrS

    filename = '../sshFiles/RTI_2D0020.dat'
    yt.set_log_level(30)

    dat1 = yt.load(fileName, unit_system='code')

    level = 2
    densityL2 = dat1.covering_grid(level, left_edge=dat1.domain_left_edge, dims = dat1.domain_dimensions * dat1.refine_by**level)['density']
    velocityL2 = dat1.covering_grid(level, left_edge=dat1.domain_left_edge, dims = dat1.domain_dimensions * dat1.refine_by**level)['velocity_y']
    yGridL2 = dat1.covering_grid(level, left_edge=dat1.domain_left_edge, dims = dat1.domain_dimensions * dat1.refine_by**level)['y']

    xnum, _, _ = np.shape(densityL2)

    rhoBubble = densityL2[0, :, 0].value
    rhoSpike = densityL2[int(xnum/2), :, 0].value
    gradRhoBubble = np.gradient(rhoBubble, edge_order=2)
    gradRhoSpike = np.gradient(rhoSpike, edge_order=2)

    indexBubbleMax = np.argmax(gradRhoBubble)
    indexSpikeMax = np.argmax(gradRhoSpike)

    UB = velocityL2[0, indexBubbleMax, 0].value
    US = velocityL2[int(xnum/2), indexSpikeMax, 0].value

    hB = yGridL2[0, indexBubbleMax, 0].value
    hS = yGridL2[int(xnum/2), indexSpikeMax, 0].value

    return hB, UB, hS, US


def getBubbleSpikeData(folderName: str, Atwoods: float, lamba: float, y0: float):
    # store bubble/spike heigts and velocities from .dat files in file bubbleSpikeData.npy. The data is store in an array in file in folderName.
    #   The columns of the array contain the bubble height, bubble velocity, spike height, spike velocity
    # Input: 
    #   folderName: folder in which all the .dat files are present
    # Outputs: 
    #   data: the same matrix that is stored is returned as well

    print(f'Retrieving bubble/spike data from {folderName}')

    filesInFolderTemp = os.listdir(folderName)
    filesInFolder = []
    for file in filesInFolderTemp:
        if '.dat' in file:
            filesInFolder.append(file)
    filesInFolder = sorted(filesInFolder)

    nbFiles = len(filesInFolder)

    data = np.zeros((nbFiles, 4))
    for index, file in enumerate(filesInFolder):
        hB, UB, hS, US = getBSHeigthVelocity(folderName + '/' + file)
        data[index, 0] = (hB - y0)/lamba
        data[index, 1] = UB/np.sqrt(Atwoods*lamba/(1.0 + Atwoods))
        data[index, 2] = abs(hS - y0)/lamba
        data[index, 3] = abs(US)/np.sqrt(Atwoods*lamba/(1.0 + Atwoods))
    
    np.save(folderName + '/' + 'bubbleSpikeData', data)

    return data

def plotBubbleSpikeData(file: str, tauMax: int = 6):
    # Plot the BubbleSpikeData 
    data = np.load(file)
    s = np.shape(data)
    tauRange = np.linspace(0, tauMax, s[0])

    plt.figure(1)

    plt.subplot(2, 2, 1)
    plt.plot(tauRange, data[:, 1], '.-')
    plt.xlabel(r'\tau')
    plt.ylabel(r'Fr_B')

    plt.subplot(2, 2, 3)
    plt.plot(tauRange, data[:, 0], '.-')
    plt.xlabel(r'\tau')
    plt.ylabel(r'h_B / lambda')

    plt.subplot(2, 2, 2)
    plt.plot(tauRange, data[:, 3], '.-')
    plt.xlabel(r'\tau')
    plt.ylabel(r'Fr_S')

    plt.subplot(2, 2, 4)
    plt.plot(tauRange, data[:, 2], '.-')
    plt.xlabel(r'\tau')
    plt.ylabel(r'h_S / lambda')

    slashIndex = file.rfind('/')
    folderName = file[0:slashIndex+1]

    plt.savefig(folderName + 'BubbleSpikeHeightVelocity.png')
    plt.show()


'''
if __name__ == '__main__':
    # The folder in which all the .dat all stored should be passed as a command line argument

    if len(sys.argv) != 5:
        raise ValueError('Command line arguments: 1. path to folder with .dat files. 2. Atwoods number. 3. lambda. 4. y0')
    Atwoods = float(sys.argv[2])
    lamba = float(sys.argv[3])
    y0 = float(sys.argv[4])
    folder = sys.argv[1]

    getBubbleSpikeData(folder, Atwoods, lamba, y0)
    plotBubbleSpikeData(folder + 'bubbleSpikeData.npy')
'''
