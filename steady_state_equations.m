% Rohan Vemu
%% Parameters for the Model
p_r = 0.5;
Lux_R = 0.02;
AHL = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2].*10^6;
delta_r = 0.0031;
alpha_TX_gfp = 0.05;
k_r = 1.35e-5;
n_1 = 1;
delta_TX_gfp = 0.2;
alpha_gfp = 2;
delta_gfp = 4;
%% Steady State Equations
R = (p_r*(Lux_R^2)*(AHL.^2))./delta_r;

TX_gfp = (alpha_TX_gfp.*(R./k_r).^n_1)./((1+(R./k_r).^n_1))./delta_TX_gfp;

GFP = alpha_gfp.*TX_gfp./delta_gfp;
%% Plot each of the steady states 
figure(1)
subplot(3,1,1)
semilogx(AHL,R,"-o","LineWidth", 2)
grid on 
grid minor
xlabel("AHL concentration (uM)")
ylabel("R Concentration (uM)")

subplot(3,1,2)
semilogx(AHL,TX_gfp,"-o","LineWidth", 2)
grid on 
grid minor
xlabel("AHL concentration (uM)")
ylabel("TX_{gfp} Concentration (uM)")

subplot(3,1,3)
semilogx(AHL,GFP,"-o","LineWidth", 2)
grid on 
grid minor
xlabel("AHL concentration (uM)")
ylabel("GFP Concentration (uM)")
