![formation of virtual cathode simulated by ESPIC_1D] (pic/fig_step_950.png)

# ESPIC_1D_VC

this is a **1D electrostatc pic** code for simulation of formation of virtual cathode in electron beam passing short circiut diode. 
for more information about pic you can see [wiki-pic](https://en.wikipedia.org/wiki/Particle-in-cell).

this code inspired from [Birdsall-Book] (https://www.taylorfrancis.com/books/mono/10.1201/9781315275048/plasma-physics-via-computer-simulation-langdon-birdsall).

# Requirements and Instalation
you need **numpy, matplotlib and tqdm** .if you haven't them simply install by `pip install numpy matplotlib tqdm`

# Running
for running code, after cloning this repo and cd to parent directory, simply run:
    `python e_beam_virtual_cathod.py`
also you can change parameters like length, time step, number of steps, plasma frequency, ... in this file and see their effect.
results by default will saved in ./results directory. 
