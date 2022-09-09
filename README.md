# ISS_Microgrids

This repo contains example simulation code for the paper 'Increasing the region of attraction in DC microgrids'. 

The simulation plots included in the paper can be replicated by running 'DCnetwork_nonIdealSim.m'. This simulation 
environment makes use of Simscape Electrical Toolbox to simulate the behaviour of switching voltage generation at the nodes.
Changing the value of Vc (line 15) between 0/24 will generate the results in the non-ideal/ideal cases, respectively. Ensure 
that Simulink and Simscape Electrical Toolbox are installed before running

An additional simulation using ideal voltage generation is included and can be run by calling 'DCnetwork_dealSim.m'.

The region of attraction plots can be reproduced by downloading 'DC_pH_1node.m' and 'RoA_1node.m' to the same directory. 
Run the section in 'RoA_1node.m' sequentially to generate each of the figures within the manuscript.
