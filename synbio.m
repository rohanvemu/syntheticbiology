%Rohan Vemu, BE310, Synthetic Biology 

function dy=synbio(t, y, i, j)

p_r = 0.5;
Lux_R = 0.02;
% Lux_R = linspace(0, 0.1, 21);
AHL = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2].*10^6;
AHL = logspace(-4,4,50);
delta_r = 0.0231;
% delta_r = linspace(0, 0.1, 201);
alpha_TX_gfp = 0.05;
k_r = 1.35e-5;
n_1 = 1;
delta_TX_gfp = 0.2;
alpha_gfp = 2;
delta_gfp = 4;

dy=zeros(3,1);

%preset constants, converted into the same units

%ODEs corresponding to R-->y(1), TX_gfp-->y(2), GFP-->y(3) 

dy(1) = p_r * Lux_R(j).^2*AHL(i).^2 - delta_r*y(1);

dy(2) = ((alpha_TX_gfp*(y(1)./k_r).^n_1) / (1 + (y(1)./k_r).^n_1)) -delta_TX_gfp*(y(2));

dy(3) = alpha_gfp*y(2) - delta_gfp*y(3);

end