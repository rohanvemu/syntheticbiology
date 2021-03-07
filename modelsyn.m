%Rohan Vemu, BE310, Synthetic Biology
%% Setting Initial Conditions and ODE Parameters
time= 120;%   in minutes
param = [0 0 0];%intial concentrations
reltol = 1e-6;
abstol = ones(1, 3) * 1e-5;
options=odeset('RelTol',1e-3,'AbsTol',abstol);
logspaceAHL = logspace(-4,4, 50);
AHL = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2].*10^6;
fluoro_store = zeros(3, 9);
%% Importing Compiled Group Data 
data = readtable('straindata.xlsx');
S1_fluoro = table2array(data(:,3));
S2_fluoro = table2array(data(:,4));
errorS1 = [0.044093821, 0.041008208, 0.051483377, 0.043257412, 0.009624476, 0.020578197, 0.027755637, 0.007059762, 0.014990832];
errorS2 = [0.023667161, 0.018169387, 0.014325166, 0.023183471, 0.007962266, 0.000575392, 0.000372541, 0.00056127, 0.000222294];
%% Modeling ODEs and Storing Final GFP Expression
for i = 1:length(logspaceAHL)
func_handle = @(t,y)synbio(t, y, i, 1);
[T,Y]=ode45(func_handle,[0 time],param,options);
fluoro_store(:, i) = Y(end, :)';
end
normalized_store = fluoro_store(3, :) / max(fluoro_store(3, :));
%% Calculating RMSE Values
% Root Mean Squared Error Formula
RMSE_S1 = 0.142;
RMSE_S2 = 0.071;
%% Plotting Fluorescence with RMSE
figure(1)
r = errorbar(logspaceAHL, normalized_store, zeros(1, length(logspaceAHL)), '--o');
set(gca,'XScale','log')
hold on 
v = errorbar(AHL,S1_fluoro,errorS1, '--o');
set(gca,'XScale','log')
set(gca,'FontSize',14)
ylim([0, 1])
set(r(1),'linewidth',1);
set(v(1),'linewidth',1.5);
legend('Model Data',['Strain 1 RMSE= ',num2str(round(RMSE_S1, 3))], 'Location', 'southeast')
xlabel("AHL Concentration [uM]")
ylabel("Normalized GFP Fluorescence")
grid on 
grid minor 
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