# Dynamin Cluster Simulation
## System Description

1. Overview
    
    This project involves developing a spatial reaction-diffusion model that recapitulates the kinetics of dynamin assembly into clusters on the membrane. The model includes the following components:
    
    - D_sol: dynamin in solution
    - R_mem: a membrane recruiter uniformly distributed on the membrane
    - A_clus: an 'activating' recruiter that is localized to the membrane cluster region (a circular region in the center of the membrane)
    - D_mem: dynamin localized to the membrane after binding the recruiter R_mem
    - D_clus: dynamin localized to the cluster site
    
    The input files for the pre-stimulation are located in the PRE folder, while the input files for the post-stimulation are located in the POST folder.
    

## Simulation Procedure

To perform the simulations, follow these steps:

1. Clone the nerdss_development repo and checkout to the dynamin_development branch.
2. Compile the nerdss executable.
3. Submit the initial simulation to reach the steady state with the dynamin molecules bound on the membrane using the command: `./nerdss -f init.inp`.
4. Then submit the restart simulation with all the reactions turned on with the restart.dat file as the input from the last step using the command: `./nerdss -a add.inp -r restart.dat`.
5. Run 48 simulations for PRE and POST systems separately, then take the average to get the dynamic.