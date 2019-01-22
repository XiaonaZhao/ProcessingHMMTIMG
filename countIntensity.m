origin
s = zeros(row, 1);
s0 = zeros(row, 1);
c = zeros(row, 1);

for n = 1:row
% mask = imread('K:\20181227_MoS2_CH18-SH\recount_IMG\mask_1227_z1_A6.tif');
mask = imread(fullfile(fileFolder1, fileNames1{n}));
mask = ~mask;
% mask1 = imread('K:\20181227_MoS2_CH18-SH\recount_IMG\mask_1227_z3_B4_bg1.tif');
% mask1 = ~mask1; % background1 mask
% mask2 = imread('K:\20181227_MoS2_CH18-SH\recount_IMG\mask_1227_z3_B4_bg2.tif');
% mask2 = ~mask2; % background2 mask
% img =  imread('K:\20181227_MoS2_CH18-SH\recount_IMG\mono_1227_z1_A6_726.tif');
img =  imread(fullfile(fileFolder2, fileNames2{n}));
img = double(img);
img0 = img.*mask;
% img1 = img.*mask1; % background1
% img2 = img.*mask2; % background2
% s0 = sum(img0(:))
s0(n) = sum(img0(:));
% bg1 = sum(img1(:));
% bg2 = sum(img2(:));
% s0 = sum(img0(:))-(bg1+bg2)/2
% c = sum(mask(:));
c(n)= sum(mask(:));

% I still haven't solve the problem of automatically shifting baseline.
% origin = 3.77313;

% s = s0 - c*origin; % shifting baseline.
% s = -s
s(n) = -(s0(n) - c(n)*origin(n)); % shifting baseline.
end
