"""
this file (e_beam_virtual_cathod) is for simulating virtual cathod formation!
you can alter input file like L,v0, NP, ...
"""
from constants import *

from generation import emiter_step

from pic_simulator import Sim


# inpyt params: ______________________________________________
PARTS_GEN = emiter_step

Omega_p = 2.0e8             # non normalized plasma frequency ~ 9 sqrt(n0)

CELLS_NUM = 100             # number of cells
#node_num = cells_num + 1       # number of node points
PARTS_NUM = 20
SIZE =  2 * (C/Omega_p)             # size of the system in meter

STEPS = 1000
PLT_STP = int(STEPS / 20) # plotting periode

dX = SIZE / CELLS_NUM       # distance between nodes (spatial step)
    
V0 = 0.2 * C                # 0.2c ~ 10keV , 0.4c~50keV
# verbose get more info of proccess
VERBOSE = False #True

# address of directory for saving plots
RESULT_DIR = "./results"
# Plasma and Beam parameter
Omega_pp= Omega_p * ( (m_e/m_p)**(0.5) )               # plasma freq for ions

dT = 1/(100*Omega_p)                # dt timestep (should be: omegaP * dt < 0.3)

# PERTURBATION
MODE = None # None for no perturbation  
AMPL = None #SIZE/5

# ____________________________________________________________

# instantiate Sim() class for simulation env:
sim = Sim(size = SIZE, cells_num = CELLS_NUM, time_step = dT,
        parts_num = PARTS_NUM,
        bound = [0,0], 
        dir = RESULT_DIR, verbose = VERBOSE)

# running PIC main loop: 
sim.run(parts_gen = PARTS_GEN, steps = STEPS, omega_p = Omega_p,
        mode = MODE , ampl = AMPL, v0 = V0,
        plot = True, plt_stp = PLT_STP)
     