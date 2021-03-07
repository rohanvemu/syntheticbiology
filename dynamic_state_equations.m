%Rohan Vemu
function [R, TXGFP, GFP, time_array] = dynamic_state_equations(AHL)
%% Parameters for the Model
p_r = 0.5;
Lux_R = 0.0001;
AHL = AHL;
delta_r = 0.0031;
alpha_TX_gfp = 0.05;
k_r = 1.35e-5;
n_1 = 1;
delta_TX_gfp = 0.2;
alpha_gfp = 2;
delta_gfp = 4;
%% Dynamic Equations 
dt = 0.1;  %define time step for Euler's method
t_end = 24*60;
time = 0:dt:t_end;

R = zeros(length(AHL),length(time)); % intialize R
TXGFP = zeros(length(AHL),length(time)); % intialize TFGFP
GFP = zeros(length(AHL),length(time)); % intialize GFP

for j=1:length(AHL)
    
for i = 1:length(time)-1
    
    dR_dt=p_r.*(Lux_R.^2)*(AHL(j).^2)-delta_r.*R(j,i);
    R(j,i+1) = R(j,i) + dt.*dR_dt;
    
    dTXGFP_dt=(alpha_TX_gfp.*(R(j,i)./k_r).^n_1)./(1+(R(j,i)/k_r).^n_1)-delta_TX_gfp*TXGFP(j,i);
    TXGFP(j,i+1)= TXGFP(j,i) + dt.*dTXGFP_dt;
    
    dGFP_dt=alpha_gfp*TXGFP(j,i) -delta_gfp*GFP(j,i);
    GFP(j,i+1)= GFP(j,i) + dt.*dGFP_dt;
    
end

end
%%
length_time = size(time, 2)-1;
span = 0:24;
time_array =(1 + (span.*length_time./24));
end

