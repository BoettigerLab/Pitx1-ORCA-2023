'''
 Polymer simulations

These simulations generate TADs without loop extrusion, (relying instead on weak sticky interactions between monomers) as shown and described in Extended Data Figure 10 and the associated text.
The four simuluations are identical except for their interactionMatrix, which is specified here: 
Swap which line is commented out (on line 53) to switch between the simulations.

interactionMatrix = np.array([[.3,.2,.2,.2,.2,.2],[.2,.3,.2,.2,.2,.2],[.2,.2,.3,.2,.2,.2],[.2,.2,.2,.3,.2,.2],[.2,.2,.2,.2,.3,.2],[.2,.2,.2,.2,.2,.8]])  # Weak TADs
# interactionMatrix = np.array([[.3,0,0,0,0,0],[0,.3,0,0,0,0],[0,0,.3,0,0,0],[0,0,0,.3,0,0],[0,0,0,0,.3,0],[0,0,0,0,0,.8]])  # Strong TADs
# interactionMatrix = np.array([[.3,0,0,0,0,0],[0,.3,0,0,0,0],[0,0,.3,0,0,0],[0,0,0,.3,0,0],[0,0,0,0,.3,0],[0,0,0,0,0,.5]])  # Weak E-clustering
# interactionMatrix = np.array([[.3,0,0,0,0,0],[0,.3,0,0,0,0],[0,0,.3,0,0,0],[0,0,0,.3,0,0],[0,0,0,0,.3,0],[0,0,0,0,0,1.1]])  # Strong E-clustering

Additional information can be found in the methods. 
'''
saveFolder = 'Updated/SavePath/'

import sys
import os
import numpy as np
import numpy.matlib
import h5py
import ast
import pandas as pd
import math

sys.path.append("C:\Shared\polychrom-shared")
sys.path.append("C:\Shared\polychrom-shared\Simulations")
from LEBondUpdater import bondUpdater
import polychrom
from polychrom.starting_conformations import grow_cubic
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
from polychrom.simulation import Simulation
from polychrom import polymerutils
from polychrom import forces
from polychrom import forcekits
import time

# general parameters 
trajectoryLength = 2000 # time duration of simulation 
density = 0.1 # chromatin density (achieved using periodic boundary conditions)
   
#  ==========Extrusion sim parameters====================
# there is probably a more elegant way to read in text values than ast.literal_eval, but this works.  
N1 = 1000 # Number of monomers in the polymer
M = 100  # Replicate polymers 
N = N1 * M # number of monomers in the full simulation 
LEFNum =0
LIFETIME = 100  # Not used, kept here to observe parallele structure and generality. 
ctcfSites = np.array([0,1000]) # Not used, kept here to observe parallele structure and generality. 
nCTCF = np.shape(ctcfSites)[0] # Not used, kept here to observe parallele structure and generality. 
ctcfDir = np.zeros(nCTCF) # 0 is bidirectional, 1 is right 2 is left.  Not used, kept here to observe parallele structure and generality. 
ctcfCapture = 0.99*np.ones(nCTCF) #  Not used, kept here to observe parallele structure and generality. 
ctcfRelease =0.01*np.ones(nCTCF)  # % Not used, kept here to observe parallele structure and generality. 
loadProb = 0  # not used, -- keep cohesin loading probability uniform
interactionMatrix = np.array([[.3,.2,.2,.2,.2,.2],[.2,.3,.2,.2,.2,.2],[.2,.2,.3,.2,.2,.2],[.2,.2,.2,.3,.2,.2],[.2,.2,.2,.2,.3,.2],[.2,.2,.2,.2,.2,.8]])  # Weak TADs
# interactionMatrix = np.array([[.3,0,0,0,0,0],[0,.3,0,0,0,0],[0,0,.3,0,0,0],[0,0,0,.3,0,0],[0,0,0,0,.3,0],[0,0,0,0,0,.8]])  # Strong TADs
# interactionMatrix = np.array([[.3,0,0,0,0,0],[0,.3,0,0,0,0],[0,0,.3,0,0,0],[0,0,0,.3,0,0],[0,0,0,0,.3,0],[0,0,0,0,0,.5]])  # Weak E-clustering
# interactionMatrix = np.array([[.3,0,0,0,0,0],[0,.3,0,0,0,0],[0,0,.3,0,0,0],[0,0,0,.3,0,0],[0,0,0,0,.3,0],[0,0,0,0,0,1.1]])  # Strong E-clustering


# specify TADs
oneChainMonomerTypes = np.zeros(N1).astype(int)
oneChainMonomerTypes[0:200] = 0
oneChainMonomerTypes[200:400] = 1
oneChainMonomerTypes[400:600] = 2
oneChainMonomerTypes[600:800]= 3
oneChainMonomerTypes[800:1000] =4
w = 5
oneChainMonomerTypes[200-w:200+w] = 5
oneChainMonomerTypes[400-w:400+w] = 5
oneChainMonomerTypes[600-w:600+w] = 5
oneChainMonomerTypes[800-w:800+w] = 5

# load prob
if loadProb == 0:
    loadProb = np.ones([1,N1])  # uniform loading probability
    loadProb = numpy.matlib.repmat(loadProb,1,M) # need to replicate and renormalize
    loadProb = loadProb/np.sum(loadProb) 

if len(oneChainMonomerTypes) != N1:
    oneChainMonomerTypes = np.zeros(N1).astype(int)

if not os.path.exists(saveFolder):
    os.mkdir(saveFolder)
lefPosFile = saveFolder + "LEFPos.h5"


# less common parameters
attraction_radius = 1.5
num_chains = M  # simulation uses some equivalent chains  (5 in a real sim)
MDstepsPerCohesinStep = 800
smcBondWiggleDist = 0.2
smcBondDist = 0.5
saveEveryBlocks = 10   # save every 10 blocks (saving every block is now too much almost)
restartSimulationEveryBlocks = 100


# =====  Skip ahead to 3D simulation 
# (this part is just kept to make it easy to add back in the loop extrusion if desired)
# ======
# 1D extrusion (not used, kept to show symmetry in the code construction and generalizability. 
import polychrom.lib.extrusion1Dv2 as ex1D # 1D classes 
ctcfLeftRelease = {}
ctcfRightRelease = {}
ctcfLeftCapture = {}
ctcfRightCapture = {}
# should modify this to allow directionality
for i in range(M): # loop over chains (this variable needs a better name Max)
    for t in range(len(ctcfSites)):
        pos = i * N1 + ctcfSites[t] 
        if ctcfDir[t] == 0:
            ctcfLeftCapture[pos] = ctcfCapture[t]  # if random [0,1] is less than this, capture
            ctcfLeftRelease[pos] = ctcfRelease[t]  # if random [0,1] is less than this, release
            ctcfRightCapture[pos] = ctcfCapture[t]
            ctcfRightRelease[pos] = ctcfRelease[t]
        elif ctcfDir[t] == 1: # stop Cohesin moving toward the right  
            ctcfLeftCapture[pos] = 0  
            ctcfLeftRelease[pos] = 1  
            ctcfRightCapture[pos] = ctcfCapture[t]
            ctcfRightRelease[pos] = ctcfRelease[t]
        elif ctcfDir[t] == 2:
            ctcfLeftCapture[pos] = ctcfCapture[t]  # if random [0,1] is less than this, capture
            ctcfLeftRelease[pos] = ctcfRelease[t]  # if random [0,1] is less than this, release
            ctcfRightCapture[pos] = 0
            ctcfRightRelease[pos] = 1
args = {}
args["ctcfRelease"] = {-1:ctcfLeftRelease, 1:ctcfRightRelease}
args["ctcfCapture"] = {-1:ctcfLeftCapture, 1:ctcfRightCapture}        
args["N"] = N 
args["LIFETIME"] = LIFETIME
args["LIFETIME_STALLED"] = LIFETIME  # no change in lifetime when stalled 
occupied = np.zeros(N)
occupied[0] = 1  
occupied[-1] = 1 
cohesins = []
for i in range(LEFNum):
    ex1D.loadOneFromDist(cohesins,occupied, args,loadProb) # load the cohesins 
with h5py.File(lefPosFile, mode='w') as myfile:
    dset = myfile.create_dataset("positions", 
                                 shape=(trajectoryLength, LEFNum, 2), 
                                 dtype=np.int32, 
                                 compression="gzip")
    steps = 100    #
    bins = np.linspace(0, trajectoryLength, steps, dtype=int) # chunks boundaries 
    for st,end in zip(bins[:-1], bins[1:]):
        cur = []
        for i in range(st, end):
            ex1D.translocate(cohesins, occupied, args,loadProb)  # actual step of LEF dynamics 
            positions = [(cohesin.left.pos, cohesin.right.pos) for cohesin in cohesins]
            cur.append(positions)  # appending current positions to an array 
        cur = np.array(cur)  # when we finished a block of positions, save it to HDF5 
        dset[st:end] = cur
    myfile.attrs["N"] = N
    myfile.attrs["LEFNum"] = LEFNum
trajectory_file = h5py.File(lefPosFile, mode='r')
LEFNum = trajectory_file.attrs["LEFNum"]  # number of LEFs
LEFpositions = trajectory_file["positions"]  # array of LEF positions  
steps = MDstepsPerCohesinStep 
Nframes = LEFpositions.shape[0] 
print(f'Length of the saved trajectory: {Nframes}')
block = 0  # starting block 
assert (Nframes % restartSimulationEveryBlocks) == 0 
assert (restartSimulationEveryBlocks % saveEveryBlocks) == 0
savesPerSim = restartSimulationEveryBlocks // saveEveryBlocks
simInitsTotal  = (Nframes) // restartSimulationEveryBlocks
if len(oneChainMonomerTypes) != N:
    monomerTypes = np.tile(oneChainMonomerTypes, num_chains)
else:
    monomerTypes = oneChainMonomerTypes    
N_chain = len(oneChainMonomerTypes)  
N = len(monomerTypes)
print(f'N_chain: {N_chain}')  # ~8000 in a real sim
print(f'N: {N}')   # ~40000 in a real sim
N_traj = trajectory_file.attrs["N"]
print(f'N_traj: {N_traj}')
assert N == trajectory_file.attrs["N"]
print(f'Nframes: {Nframes}')
print(f'simInitsTotal: {simInitsTotal}')


#==============================================================#
#                  RUN 3D simulation                              #
#==============================================================#
milker = bondUpdater(LEFpositions)
data = grow_cubic(N,int((N/(density*1.2))**0.333))  # starting conformation
reporter = HDF5Reporter(folder=saveFolder, max_data_length=50)
PBC_width = (N/density)**0.333
chains = [(N_chain*(k),N_chain*(k+1),0) for k in range(num_chains)]

for iteration in range(simInitsTotal):
    a = Simulation(N=N, 
                   error_tol=0.01, 
                   collision_rate=0.01, 
                   integrator ="variableLangevin", 
                   platform="cuda",
                   GPU = "0", 
                   PBCbox=(PBC_width, PBC_width, PBC_width),
                   reporters=[reporter],
                   precision="mixed")  # platform="CPU", 
    a.set_data(data)
    a.add_force(
        polychrom.forcekits.polymer_chains(
            a,
            chains=chains,
            nonbonded_force_func=polychrom.forces.heteropolymer_SSW,
            nonbonded_force_kwargs={
                'attractionEnergy': .05,  # base attraction energy for all monomers
                'attractionRadius': attraction_radius,
                'interactionMatrix': interactionMatrix,
                'monomerTypes': monomerTypes,
                'extraHardParticlesIdxs': []
            },
            bond_force_kwargs={
                'bondLength': 1,
                'bondWiggleDistance': 0.05
            },
            angle_force_kwargs={
                'k': 0.5# 1.5
            }
        )
    )
    
    # ------------ initializing milker; adding bonds ---------
    # copied from addBond
    kbond = a.kbondScalingFactor / (smcBondWiggleDist ** 2)
    bondDist = smcBondDist * a.length_scale

    activeParams = {"length":bondDist,"k":kbond}
    inactiveParams = {"length":bondDist, "k":0}
    milker.setParams(activeParams, inactiveParams)
     
    # this step actually puts all bonds in and sets first bonds to be what they should be
    milker.setup(bondForce=a.force_dict['harmonic_bonds'],
                blocks=restartSimulationEveryBlocks)

    # If your simulation does not start, consider using energy minimization below
    if iteration == 0:
        a.local_energy_minimization() 
    else:
        a._apply_forces()
    
    for i in range(restartSimulationEveryBlocks):        
        if i % saveEveryBlocks == (saveEveryBlocks - 1):  
            a.do_block(steps=steps)
        else:
            a.integrator.step(steps)  # do steps without getting the positions from the GPU (faster)
        if i < restartSimulationEveryBlocks - 1: 
            curBonds, pastBonds = milker.step(a.context)  # this updates bonds. You can do something with bonds here
    data = a.get_data()  # save data and step, and delete the simulation
    del a
    
    reporter.blocks_only = True  # Write output hdf5-files only for blocks
    
    time.sleep(0.2)  # wait 200ms for sanity (to let garbage collector do its magic)

reporter.dump_data()