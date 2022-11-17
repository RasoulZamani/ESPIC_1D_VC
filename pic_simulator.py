# import libs: ________________________________________________________
import os
import numpy as np
from constants import *
import matplotlib.pyplot as plt
from tqdm import tqdm 
from generation import emiter_step
from cycle import density, sor_solver, field_on_nodes, \
field_on_particles, rewind, move_particles

# _____________________________________________________________________
class Sim():
    """class for simulation of pic
    """
    def __init__(self, size: float, cells_num, time_step, parts_num, bound, dir, verbose = False):
        self.size = size
        self.cells_num = cells_num
        self.parts_num = parts_num  # number of particles emits in dt / or total particles 
        self.time_step = time_step
        self.bound     = bound      # boundary condition: tupple or list (phi_0,phi_L)
        self.dx = size / cells_num
        self.dir = dir
        self.verbose = verbose
        
        # creating directory for savuing results
        dir_path = os.path.join(self.dir)
        os.makedirs(dir_path, exist_ok=True)
        # clear possible files from previouse run:
        files = os.listdir(dir_path) # list of existing file in dir
        if len(files)>0:
            for f in files:
                os.remove(os.path.join(dir_path,f))
       
        if self.verbose:
            print(f"\n directory:  {dir}  created for saving result\n")
    

    def run(self,
            parts_gen, # particle generator such two stream or electron beam
            steps,
            omega_p,
            mode , ampl ,
            v0 ,
            plot, plt_stp ):
        
        
        if self.verbose:
            # ..... TODO ... : print other parameters too
            print( " particles in t=0 generated with ...")
            print(f"\nnumber of cells is: {self.cells_num}, length of env: {self.size}m",
                f"\nlength of each grid(dX): {self.dx}, time step: {self.time_step}, number of steps:{steps}",
                f"\nplasma freq for e:{omega_p}, beam velocity: {v0} m/s")
        
        self.particles = []
        # main loop of PIC ___________________________________________
        for step in tqdm(range(steps), desc="pic steps"):
            if self.verbose:
                print(f"\n************   we are in step {step} ********************** \n")
            
            # generating particle with v0 and Wp in forst node:
            self.emits = parts_gen(
                        omega_p = omega_p, 
                        Ne = self.parts_num,
                        #mode = mode , ampl = ampl,
                        dx = self.dx, 
                        v0 = v0)
            # adding new emited particles in this time step:
            self.particles.extend(self.emits)
            if self.verbose: 
                print("\nnew emits particles was added")
            
            rho = density(self.particles, self.cells_num, self.dx)
            if self.verbose: 
                print("\nwheiting dencity on grids from particle position (rho) by dencity funcion done! \n")
    
            PHI = sor_solver(rho, self.cells_num, self.dx, self.bound[0], self.bound[1])
            if self.verbose:
                print(f"\nPoisson solved by SOR method")

                
            EFIELDn = field_on_nodes(PHI, self.cells_num, self.dx)
            if self.verbose:
                print("\ncalculating field from potential on grids (efild) done!\n")

                
            EFIELD = field_on_particles(EFIELDn, self.particles, self.cells_num, self.dx)
            if self.verbose:       
                print("\ncalculating E field on particles done! \n")

            if step == 0:
                self.particles = rewind(-1, EFIELD, self.particles, dt=self.time_step)

            self.particles = move_particles(EFIELD, self.particles, self.size, dt= self.time_step)
            if self.verbose:
                print("\nmoving of particles done! \n")

            # write to file
            #  ..... TODO ....
    
            # plotting _______________________________________________

            if (plot) & (step % plt_stp == 0):
                n_parts = len(self.particles)
                xi = []
                vi = []
                for i in range(n_parts):
                    if self.particles[i].mv:
                        xi.append(self.particles[i].x)
                        vi.append(self.particles[i].v)
        
                fig = plt.figure()      
                plt.scatter(x=xi, y=vi)
                plt.xlabel(" x [position] ")
                plt.ylabel(" v (velocity) ")
                plt.title(f" phase space (x-v) of particles in step {step} ")
                fig.savefig(r"results/fig_"+f"step_{step}")
        
        
        
                # ...... TODO ....... 
                # make plot more pretty later!
