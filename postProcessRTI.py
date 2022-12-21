import numpy as np
import matplotlib.pyplot as plt
import yt
import os
import sys
from typing import Final, Literal
from scipy.constants import m_p, k, R, N_A

# plt.rcParams['text.usetex'] = True # Use latex formatting in matplotlib figures 

def getParfileValues(Atwoods: float, Reynolds: float, lamba: float = 0.4, rholight: float = 1.0, taumax: float = 6, he_abundance = 0.1, g = 1.0):
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

    m_pCGS: Final[float] = m_p*1e3
    kCGS: Final[float] = k*1e7
    RCGS: Final[float] = R*1e7

    # Dimensions
    unit_time: float = np.sqrt(lamba/Atwoods) # which lambda?
    unit_velocity = np.sqrt(lamba*Atwoods/(1.0+Atwoods))
    unit_numberdensity: float =1e20
    unit_density: float = (1 + 4*he_abundance)*m_pCGS*unit_numberdensity
    unit_length = unit_velocity*unit_time
    unit_mass = unit_density*np.power(unit_length, 3)
    unit_pressure = unit_density*np.power(unit_velocity, 2)
    unit_temperature = unit_pressure/((2 + 3*he_abundance)*unit_numberdensity*kCGS)

    # Parameters with dimension
    gCGS = g*unit_length/(np.power(unit_time, 2))  # in cm/s^2
    rholightCGS: float = rholight*unit_density
    rhodenseCGS: float =rholightCGS*(1.0+Atwoods)/(1.0-Atwoods)  # in g/cm^3
    lambaCGS = lamba*unit_length  # in cm
    vc_muCGS = lambaCGS*np.sqrt(Atwoods*lambaCGS*gCGS/(Atwoods+1.0))/(2*Reynolds/(rholightCGS + rhodenseCGS))  # in g/(cm s)
    time_maxCGS = taumax/np.sqrt(Atwoods*g/(lamba))  # in s
    heatCapCGS = (5/2)*RCGS*unit_numberdensity/N_A  # in erg/K

    # Dimensionless variables
    vc_mu = vc_muCGS*unit_length*unit_time/unit_mass
    heatCap = heatCapCGS*unit_temperature/(unit_pressure)
    tc_k_para = vc_mu*heatCap

    print(f'{unit_time=}')
    print(f'{unit_velocity=}')
    print(f'{unit_density=}')
    print(f'{unit_velocity=}')
    print(f'{unit_length=}')
    print(f'{unit_temperature=}')
    print(f'{vc_mu=}')
    print('------- Par file & mod_usr.t -------')
    print(f'{unit_numberdensity=}')
    print(f'{vc_muCGS=}')
    print(f'{tc_k_para=}')
    print(f'{time_maxCGS=}')

    return vc_muCGS, tc_k_para, time_maxCGS

def getBSHeigthVelocity(fileName: str, level: int = 2):
    # Get the Bubble and spike height and velocity
    # Inputs:
    #   fileName: str: name of .dat file
    #   level: int: level at which to interpolate
    # Outputs:
    #   Buuble/spike heights/velocities: double: hB, FrB, hS, FrS

    filename: Literal['../sshFiles/RTI_2D0020.dat'] = '../sshFiles/RTI_2D0020.dat'
    yt.set_log_level(30)

    dat1 = yt.load(fileName, unit_system='code')

    level: int = 2
    densityL2 = dat1.covering_grid(level, left_edge=dat1.domain_left_edge, dims = dat1.domain_dimensions * dat1.refine_by**level)['density']
    velocityL2 = dat1.covering_grid(level, left_edge=dat1.domain_left_edge, dims = dat1.domain_dimensions * dat1.refine_by**level)['velocity_y']
    yGridL2 = dat1.covering_grid(level, left_edge=dat1.domain_left_edge, dims = dat1.domain_dimensions * dat1.refine_by**level)['y']

    xnum, _, _ = np.shape(densityL2)

    rhoBubble = densityL2[0, :, 0].value
    rhoSpike = densityL2[int(xnum/2), :, 0].value
    gradRhoBubble = np.gradient(rhoBubble, edge_order=2)
    gradRhoSpike = np.gradient(rhoSpike, edge_order=2)

    indexBubbleMax: signedinteger = np.argmax(gradRhoBubble)
    indexSpikeMax: signedinteger = np.argmax(gradRhoSpike)

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
        print(f'Done with {file}')
        hB, UB, hS, US = getBSHeigthVelocity(folderName + '/' + file)
        data[index, 0] = (hB - y0)/lamba
        data[index, 1] = UB/np.sqrt(Atwoods*lamba/(1.0 + Atwoods))
        data[index, 2] = abs(hS - y0)/lamba
        data[index, 3] = abs(US)/np.sqrt(Atwoods*lamba/(1.0 + Atwoods))
    
    np.save(folderName + '/' + 'bubbleSpikeData', data)

    return data

def plotBubbleSpikeData(file: str, tauMax: int = 6, tauCutOff = None) -> None:
    # Plot the BubbleSpikeData 
    data = np.load(file)
    s = np.shape(data)
    tauRange = np.linspace(0, tauMax, s[0])

    if tauCutOff is not None:
        maxIndex = int(s[0]*tauCutOff)
        tauRange = tauRange[0:maxIndex]
        data = data[0:maxIndex, :]

    plt.figure(1, figsize=(10, 6))

    plt.subplot(2, 2, 1)
    plt.plot(tauRange, data[:, 1], '.-')
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$Fr_B$')

    plt.subplot(2, 2, 3)
    plt.plot(tauRange, data[:, 0], '.-')
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$h_B / \lambda$')

    plt.subplot(2, 2, 2)
    plt.plot(tauRange, data[:, 3], '.-')
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$Fr_S$')

    plt.subplot(2, 2, 4)
    plt.plot(tauRange, data[:, 2], '.-')
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$h_S / \lambda$')

    slashIndex: int = file.rfind('/')
    folderName: str = file[0:slashIndex+1]

    plt.savefig(folderName + 'BubbleSpikeHeightVelocity.png')
    plt.show()

def getVorticity(dataFile, imageFile, A=0.04, g=1, lamba=0.4) -> None:
    # Make a vorticity plot
    # Inputs:
    #   dataFile: .dat file of which to make a vorticity plot
    #   imageFile: file to save image in
    #   A: Atwood number
    #   g: gravity constant
    #   lamba: perturbation wavelength

    dat1 = yt.load(dataFile, unit_system='code')

    level: Final[int] = 1
    velocityXGrid = dat1.covering_grid(level, left_edge=dat1.domain_left_edge, dims = dat1.domain_dimensions * dat1.refine_by**level)['velocity_x']
    velocityYGrid = dat1.covering_grid(level, left_edge=dat1.domain_left_edge, dims = dat1.domain_dimensions * dat1.refine_by**level)['velocity_y']

    velocityX = velocityXGrid[:, :, 0].value
    velocityY = velocityYGrid[:, :, 0].value

    del velocityYGrid
    del velocityXGrid

    t1 = np.gradient(velocityX, edge_order=2, axis=1)
    t2 = np.gradient(velocityY, edge_order=2, axis=0)
    vorticity = (t1-t2)/np.sqrt(A*g/lamba)

    plt.matshow(np.rot90(vorticity))
    plt.savefig(imageFile)
    plt.show()


def plotBubbleSpikeDataComparison(file: str, file2: str, tauMax: int = 6, tauCutOff = None) -> None:
    # Plot the BubbleSpikeData 
    data1 = np.load(file)
    data2: Any = np.load(file2)
    s = np.shape(data1)
    tauRange = np.linspace(0, tauMax, s[0])

    if tauCutOff is not None:
        maxIndex = int(s[0]*tauCutOff)
        tauRange = tauRange[0:maxIndex]
        data1 = data1[0:maxIndex, :]
        data2 = data2[0:maxIndex, :]

    plt.figure(1, figsize=(10, 6))

    plt.subplot(2, 2, 1)
    plt.plot(tauRange, data1[:, 1], '.-')
    plt.plot(tauRange, data2[:, 1], '.-')
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$Fr_B$')

    plt.subplot(2, 2, 3)
    plt.plot(tauRange, data1[:, 0], '.-')
    plt.plot(tauRange, data2[:, 0], '.-')
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$h_B / \lambda$')

    plt.subplot(2, 2, 2)
    plt.plot(tauRange, data1[:, 3], '.-')
    plt.plot(tauRange, data2[:, 3], '.-')
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$Fr_S$')

    plt.subplot(2, 2, 4)
    plt.plot(tauRange, data1[:, 2], '.-')
    plt.plot(tauRange, data2[:, 2], '.-')
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$h_S / \lambda$')

    slashIndex: int = file.rfind('/')
    folderName: str = file[0:slashIndex+1]

    plt.savefig(folderName + 'BubbleSpikeHeightVelocityComparison.png')
    plt.show()



if __name__ == '__main__':
    # -- Run this to get the parameters for the parfile and mod_usr.t file
    getParfileValues(0.04, 100)


    # -- This is for plotting the bubble/spike height/velocities
    if len(sys.argv) != 5:
        raise ValueError('Command line arguments: 1. path to folder with .dat files. 2. Atwoods number. 3. lambda. 4. y0')
    Atwoods = float(sys.argv[2])
    lamba = float(sys.argv[3])
    y0 = float(sys.argv[4])
    folder = sys.argv[1]

    getBubbleSpikeData(folder, Atwoods, lamba, y0)
    plotBubbleSpikeData(folder + 'bubbleSpikeData.npy')

