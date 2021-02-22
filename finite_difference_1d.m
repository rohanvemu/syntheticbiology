%Rohan Vemu, BE310, Synthetic Biology Problem Set
%Modeling 1D Fick's Second Law via Finite Difference Approximation 
%% Defining Parameters
L = 100; %in mm
T = 60; %in minutes
D = 1; % in mm^2/min
IC_reservoir = 100; % in uM
%% Defining the Geometry
n = 101; % number of grid points
dx = L/(n-1); %setting a dx
distance = [0:dx:L]; %distance vector incremented by dx
%% Discretization of Time
stability_factor = 0.5; % must be <= 0.5 for FTCS to converge to correct solution
dt = stability_factor*(dx^2)/D; % dt that fulfills the stability criterion
time = [0:dt:T]; % time vector incremented by dt
r = D*dt/(dx^2); % r term in the Finite Difference Equation
%% Specifying Initial/Boundary Conditions
diff_store = zeros(length(time), n); %matrix to store evolving concentrations
diff_store(:, 1) = IC_reservoir; %setting the initial conditions
diff_store(:, end) = 0; %setting the BC
%% Solving the Diffusion Equation 1D
for i = 1:(length(time)-1) %need to exclude first/last element since j-1, j+1, does not exist for those
    for j = 2:(n-1)
        %perform finite difference approximation
        diff_store(i+1, j) = diff_store(i, j) + r*(diff_store(i, j+1) + diff_store(i, j-1) - 2*diff_store(i, j));
    end   
end
%% Plotting the Distance/Time Evolution of the System
%plots of concentration versus distance
thirty_min_ind = (length(time) + 1) / 2; %find index corresponding to 30 min
figure(1)
hold on 
plot(distance, diff_store(1, :))
plot(distance, diff_store(thirty_min_ind, :))
plot(distance, diff_store(end, :))
%standard plotting commands and labeling
xlabel("Distance (mm)")
ylabel("Concentration (uM)")
title("Time Evolution of Concentration")
legend('t = 0 min', 't = 30 min', 't = 60 min')

%plots of concentration versus time
figure(2)
hold on 
plot(time, diff_store(:, 0+1)) %array is constructed such that time and index are 1:1
plot(time, diff_store(:, 5+1))
plot(time, diff_store(:, 15+1))
plot(time, diff_store(:, 100+1))
%standard plotting commands and labeling
xlabel("Time (min)")
ylabel("Concentration (uM)")
title("Evolution of Concentration at a Given Distance")
legend('x = 0 mm', 'x = 5 mm', 'x = 15 mm', 'x = 100 mm')