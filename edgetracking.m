%% Importing Images
files = dir('*2.0s_Exposure*.png');
image_store=zeros(24, size(imread(files(1).name), 1), size(imread(files(1).name), 2));

for i = 1:length(files)-1
    A = rgb2gray(imread(files(i).name));
    B = rgb2gray(imread(files(end).name));
    new = A-B;
    image_store(i, :, :) = new./max(new, [], 'all');
end
imshow(files(end-1).name);
hp = impixelinfo;
[x, y] = ginput(2); %select center point and then any point on radius of plate
%% Using Line Profile to GFP Edge Distance
midpoint_j = round(y(1));
midpoint_i = round(x(1));
threshold= 0.1;
index_store=zeros(1, length(image_store(:,1,1)));

for j=1:length(image_store(:,1,1))

for k=midpoint_i:length(image_store(1,1,:))
    val=image_store(j, midpoint_j, k);
    if val<threshold
        index_store(j)=k;
        break
    end
end

end
%% Determine Pixel to Distance Ratio and Convert Values
radius_plate = 4.25;
pix_length = norm([x(2), y(2)] - [x(1), y(1)]);
pix_to_cm = radius_plate / pix_length;

pix_from_center = index_store - round(x(1));
dist_from_center = pix_from_center .* pix_to_cm;