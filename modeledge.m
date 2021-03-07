%Rohan Vemu
%% Find GFP Edge Distance for Model 
length_time = size(AHL, 3);
span = 1:24;
time_index_store = round(span.*length_time./24); %find indices in AHL/time array corresponding to each hour
store_AHL = AHL(:, :, (time_index_store)); %AHL values per hour, need to convert to GFP values
line_AHL = store_AHL(101, 101:end, :);
%% Use AHL Values and Dynamic Equations to Find GFP at Each Time
GFP_store = zeros(size(line_AHL, 2), size(line_AHL, 3));
for i = 1:size(line_AHL, 3)
[R, TXGFP, GFP, time_array] = dynamic_state_equations(line_AHL(:, :, i));
GFP_store(:, i) = GFP(:, time_array(i))/max(GFP(:, time_array(i)));
end
%% Performing Edge Detection with GFP Values
%working columnwise to do edge detection beneath a certain threshold 
threshold = 0.1;
model_index_store = zeros(1, 24);
for p = 1:size(GFP_store, 2)
for z = 1:size(GFP_store, 1)
    val=GFP_store(z, p);
    if val<threshold
        model_index_store(p)=z;
        break
    end
end
end 
%% Converting to Distance from Pixel Indices
radius_plate = 4.5;
pixel_radius = 101;
cmperpixel = radius_plate/ pixel_radius;
model_edges = model_index_store.*cmperpixel;
dist_from_center(2:end) = dist_from_center(2:end);
dist_3_center = dist_from_center(1:3:end);

NRMSE1 = RMSE(model_edges(:, 1:21), dist_from_center(:, 1:21)) / (max(dist_from_center(:, 1:21)) - min(dist_from_center(:, 1:21)))
NRMSE2 = RMSE(model_edges(:, 1:3:24), dist_3_center) / (max(dist_3_center)) - min(dist_3_center);
%%
figure(2)
r = plot(model_edges(1:end), '--o'); 
set(r(1),'linewidth',1);
hold on
r= plot(1:3:24, dist_3_center, '--o');
set(r(1),'linewidth',1);
hold on
r = plot(AHL_dist_store_s2./10, '--o');
set(r(1),'linewidth',1);
r = plot(AHL_dist_store./10, '--o');
set(r(1),'linewidth',1);
xlabel("Time (hrs)")
ylabel("Edge Distance (cm)")
grid on 
grid minor
legend("GFP Diffusion Model", ['Strain 2, NRMSE= ',num2str(round(NRMSE2, 3))], "Switchpoint AHL Front", "Total AHL Front", 'Location', 'southeast')
set(gca,'FontSize',14)


