# syntheticbiology

Code for Synthetic Biology lab module and accompanying work is presented. Peruse as you wish!

Included in this repository:
1. synbio.m --- sets up ODE functions governing genetic circuit (set up for use by ODE45)
2. modelsyn.m --- use synbio.m to model varying condition of circuit (can use ODE45 or dynamic state equations)
3. RMSE.m --- function calculating RMSE between two arrays
4. finite_difference_1d.m --- models diffusion in 1-D using finite difference
5. finite_difference_2d.m --- models diffusion in 2-D using finite difference (adjust parameters 
6. synthbio.ipynb --- normalizes and organizes pooled class data into strain_data.xlsx
7. strain_data.xlsx --- exported data from synthbio.ipynb containing model to which synbio is being fit
8. steady_state_equations.m --- steady state expressions solved for, same as taking last time point of synbio w/ ODE45 
9. dynamic_state_equations.m --- solves ODE's from synbio using Euler's method (First Order Runge Kutta)
10. modeledge.m --- runs edge tracking algorithm on 2-D diffusion data
11. edgetracking.m --- runs dge tracking algorithm on GFP plate images 


