%Rohan Vemu, Synthetic Biology, BE310
function out = RMSE(a , b)

out = sqrt(sum((a- b).^2)./length(a)); 

end