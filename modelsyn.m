%Rohan Vemu, BE310, Synthetic Biology
%% Setting Initial Conditions and ODE Parameters
time= 120;%   in minutes
param = [0 0 0];%intial concentrations
reltol = 1e-3;
abstol = ones(1, 3) * 1e-2;
options=odeset('RelTol',1e-3,'AbsTol',abstol);
AHL = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2].*10^6;
fluoro_store = zeros(3, 9);
%% Importing Compiled Group Data 
data = readtable('straindata.xlsx');
S1_fluoro = table2array(data(:,3));
S2_fluoro = table2array(data(:,4));
%% Modeling ODEs and Storing Final GFP Expression
for i = 1:length(AHL)
func_handle = @(t,y)synbio(t, y, i, 1);
[T,Y]=ode45(func_handle,[0 time],param,options);
fluoro_store(:, i) = Y(end, :)';
end
normalized_store = fluoro_store(3, :) / max(fluoro_store(3, :));
%% Calculating RMSE Values
% Root Mean Squared Error Formula
RMSE_S1 = RMSE(normalized_store, S1_fluoro');
RMSE_S2 = RMSE(normalized_store, S2_fluoro');
%% Plotting Fluorescence with RMSE
figure(1)
h = semilogx(AHL, S2_fluoro, '--o', AHL, normalized_store, '--o')
set(gca,'FontSize',14)
set(h(1),'linewidth',1.5);
set(h(2),'linewidth',1.5);
legend(['Strain 2 RMSE= ',num2str(round(RMSE_S2, 3))], 'Model Data', 'Location', 'northwest')
xlabel("AHL Concentration [M]")
ylabel("Normalized GFP Fluorescence")
grid on 
grid minor 
%% Determining How Each Constant Affects the Plot

%% Finding the Ideal Constants to Minimzie RMSE Across Both Strains
strain_number = 2;
LUX_R = linspace(0, 0.1, 21);
RMSE_store = zeros(1, length(LUX_R));
GFP_store = zeros(1, length(AHL));
for j = 1:length(LUX_R)
    for i = 1:length(AHL)
    func_handle = @(t,y)synbio(t, y, i, j);
    [T,Y]=ode45(func_handle,[0 time],param,options);

    GFP_store(i) = Y(end, 3)';
    end
    
normalized_GFP = GFP_store./max(GFP_store);    
RMSE_store(1, j) = RMSE(normalized_GFP, S1_fluoro');    
RMSE_store(2, j) = RMSE(normalized_GFP, S2_fluoro');
end
[M, I] = min(RMSE_store, [], 2)