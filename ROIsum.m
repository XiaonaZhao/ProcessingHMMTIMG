function m = ROIsum(img, mask)
% get the sum of all non-zero element in img

img = img.*mask;
m = sum(img(:));