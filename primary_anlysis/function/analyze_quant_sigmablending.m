function [sum_area,sum_stack] = analyze_quant_sigmablending(stack)
%CHECK Summary of this function goes here
%   Detailed explanation goes here
percentile_list = 75:5:95;
nsli = size(stack,3);
pixel_sli = reshape(stack,[],nsli);
%%
sli_prctile = NaN(numel(percentile_list),nsli);

for sli = 1:nsli
    pixels = pixel_sli(:,sli);
    pixels = pixels(~isnan(pixels));
    sli_prctile(:,sli) = prctile(pixels,percentile_list,"all");
end
%%
gauss_sigma = 1:1:5;
%%
sum_stack = zeros(size(stack));
for sli = 1:nsli
    hold_slice = stack(:,:,sli);
    sum_slice = zeros(size(hold_slice));
    for sigma = 1:numel(gauss_sigma)
        sum_slice = sum_slice + imgaussfilt(hold_slice,sigma);
    end
    sum_stack(:,:,sli) = sum_slice;
end
%%
sumpixel_sli = reshape(sum_stack,[],nsli);

for sli = 1:nsli
    pixels = pixel_sli(:,sli);
    pixels = sumpixel_sli(~isnan(pixels));
    sum_stack(:,:,sli) = sum_stack(:,:,sli) > prctile(pixels,75,"all");
end

sum_area = squeeze(sum(sum_stack,[1,2]));
end

