%Rohan Vemu, BE310, Synthetic Biology Problem Set
%Modeling 2D Fick's Second Law via Finite Difference Approximation 
%% Defining Parameters
radius_plate = 42.5; %in mm
radius_disk = 2.5; %in mm 
T = 24*60; %in minutes
D = 0.075; %mm^2/min
d_AHL = 4.8135e-4; %min-1
sourceconc = 10; %in uM
%% Defining the X,Y Mesh 
dx = radius_plate/100; %note dx=dy=delta p, for 201 points 
xgrid = [-radius_plate:dx:radius_plate]; % centers the grid in x around zero
ygrid = [-radius_plate:dx:radius_plate]; % centers the grid in y around zero
[X,Y] = meshgrid(xgrid, ygrid); % creates the matrix that describes the mesh using meshgrid 
%% Discretization of Time
stability_factor = 0.1; % must be <= 0.25 for FTCS to converge to correct solution
dt = stability_factor*(dx^2)/D; % dt that fulfills the stability criterion
time = [0:dt:T]; % time vector incremented by dt
r = D*dt/(dx^2); % r term in the Finite Difference Equation
%% Setting Up Initial Conditions and QC Check
AHL_initial = zeros(length(xgrid), length(ygrid));% Set initial AHL concentration to zero across the whole plate decribed by the matrix AHL_Initial
[disk_indices_row, disk_indices_col] = find(sqrt(X.^2+Y.^2)<=radius_disk); % the loop below sets the initial concentration of the disk
for i = 1:length(disk_indices_row)
 AHL_initial(disk_indices_row(i), disk_indices_col(i)) = sourceconc; %set everything in circular disk to sourceconc
end
mesh(X,Y,AHL_initial)% Draws the mesh to visualize initial conditions as a quality control step. 
title("QC Plot for Initialization")
xlabel("x distance about center of plate (mm)", 'Rotation',20)
ylabel("y distance about center of plate (mm)", 'Rotation',-30)
zlabel("Concentration (uM)")
%% Solving the Diffusion Equation in 2D
AHL = zeros(size(AHL_initial, 1), size(AHL_initial, 2), size(time, 2)); %initialization 
AHL(:, :, 1) = AHL_initial; %initialize first layer of tensor to 

for k = 1:(length(time)-1) %need to exclude last time pt
    for j = 2:(size(AHL_initial, 2) - 1) %to perform FTCS finite difference need to exclude first row/col
        for i = 2:(size(AHL_initial, 2) - 1)
            %perform finite difference approximation 
            AHL(i, j, k+1) = AHL(i, j, k) + r*(AHL(i-1, j, k)+AHL(i+1, j, k)+AHL(i, j-1, k)+AHL(i, j+1, k)-4.*AHL(i, j, k)) - dt.*AHL(i, j, k)*d_AHL;
        end
    end  
end
%% Plotting Time/Distance Evolution of the System
%2D concentration over 0, 8, 16, 24 hrs
ind_8 = floor((length(time)+1)/3); %find index corresponding to 8 hrs
ind_16 = floor(2*(length(time)+1)/3); %find index corresponding to 16 hrs

for z = [1, ind_8, ind_16, length(time)]
figure()
imagesc(xgrid, ygrid, AHL(:, :, z)) %plot the concentration map
title(['2D Concentration Profile at t= ',num2str(round(time(z)./60), 2),' hr'])
xlabel("x distance about center of plate (mm)")
ylabel("y distance about center of plate (mm)")
c = colorbar;
c.Label.String = 'Concentration (uM)'; 
%set colors/labels for the concentration maps
end

%concentration at p=10, 25mm over the 24hrs
ind_10 = floor(10/dx); %find index corresponding to 10mm 
ind_25 = floor(25/dx); %find index corresponding to 25mm 
figure()
for z = [ind_10, ind_25]
hold on 
plot(time, squeeze(AHL(101, 101+z, :)))
%center of matrix is (101, 101), moving straight up ind_10 units will get
%you to a point on the 10mm circle. expect all points on circle to have the
%same value due to random walk model and isotropic diffusion.
end
xlabel("Time (hr)")
ylabel("Concentration (um)")
legend("p=10mm", "p=25mm")
title('Concentration at a Given Radius across Time')

%want line section profile through the plate center (constant y), varying t
figure()
for z = [1, ind_8, ind_16, length(time)]
hold on 
plot(xgrid, squeeze(AHL(:, 101, z))/max(AHL(:, 101, z))) %normalize to the max value so everything can be plotted
end

title("Line Profile through Plate Center")
xlabel("x distance about center of plate (mm)")
ylabel("Normalized Concentration")
legend("t=0hr", "t=8hr", "t=16hr", "t=24hr")