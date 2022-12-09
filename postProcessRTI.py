import numpy as np
import matplotlib.pyplot as plt
import yt
import os

# plt.rcParams['text.usetex'] = True # Use latex formatting in matplotlib figures 

def getParfileValues(Atwoods, Reynolds, lamba = 0.4, rholight=1.0, taumax = 6):
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
    R = 8.31446261815324
    tc_k_para = vc_mu*5/2*R
    return vc_mu, time_max, tc_k_para


def getBSHeigthVelocity(fileName, level=2):
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

    FrB = velocityL2[0, indexBubbleMax, 0].value
    FrS = velocityL2[int(xnum/2), indexSpikeMax, 0].value

    hB = yGridL2[0, indexBubbleMax, 0].value
    hS = yGridL2[int(xnum/2), indexSpikeMax, 0].value

    return hB, FrB, hS, FrS


def getBubbleSpikeData(folderName):
    # store bubble/spike heigts and velocities from .dat files in file bubbleSpikeData.npy. The data is store in an array in file in folderName.
    #   The columns of the array contain the bubble height, bubble velocity, spike height, spike velocity
    # Input: 
    #   folderName: folder in which all the .dat files are present
    # Outputs: 
    #   data: the same matrix that is stored is returned as well

    filesInFolderTemp = os.listdir(folderName)
    filesInFolder = []
    for file in filesInFolderTemp:
        if '.dat' in file:
            filesInFolder.append(file)
    filesInFolder = sorted(filesInFolder)

    nbFiles = len(filesInFolder)

    data = np.zeros((nbFiles, 4))
    for index, file in enumerate(filesInFolder):
        hB, FrB, hS, FrS = getBSHeigthVelocity(folderName + '/' + file)
        data[index, 0] = hB
        data[index, 1] = FrB
        data[index, 2] = hS
        data[index, 3] = FrS
    
    np.save(folderName + '/' + 'bubbleSpikeData', data)

    return data

def plotBubbleSpikeData(file, tauMax):
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
    plt.ylabel(r'h_B')

    plt.subplot(2, 2, 2)
    plt.plot(tauRange, data[:, 3], '.-')
    plt.xlabel(r'\tau')
    plt.ylabel(r'Fr_S')

    plt.subplot(2, 2, 4)
    plt.plot(tauRange, data[:, 2], '.-')
    plt.xlabel(r'\tau')
    plt.ylabel(r'h_S')

    plt.show()


if __name__ == '__main__':
    getBubbleSpikeData('../sshFiles')
    plotBubbleSpikeData('../sshFiles/bubbleSpikeData.npy', 6)

