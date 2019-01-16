% main function for ProcessingHMMTIMG
clc
clear

% ATTENTION:

% There is still two problems in this programme. one is that the CURRENT
% calculated by the iLaplace doesn't match with the result from the
% electrochemical workstation. The other is that the present background
% cutting doesn't correct the bias from the charge and discharge of gold
% plate.

% Author: nonazhao@mail.ustc.edu.cn;
% Created: 14 January 2018

%%
% 1. Read .tiff files names
fileFolder = fullfile('E:\20181227_MoS2_CH18-SH\A6_zone1_10mMRu250mMPBNa_0 -0-4V_0-1VpS_2c_MoS2_uncorrosion_CH18-S-Au_HMMT100fps');
dirOutput = dir(fullfile(fileFolder,'*.tif'));
dirFile = sortObj(dirOutput);
fileNames = {dirFile.name}';
row = size(fileNames, 1);

% 2. Read one .tiff file
% 3. Use a fixed mask on the image into new one
% 4. Sum the the new image to return x
% 5. Release the original image
% 6. loop 2 to 5
intensity = zeros(row, 1);
Mask =~ imread('E:\20181227_MoS2_CH18-SH\recount_IMG\mask_1227_z1_A6.tif');
tif0 = double(imread(fullfile(fileFolder, fileNames{1})));
for n = 1: row
    tif1 = double(imread(fullfile(fileFolder, fileNames{n}))) - tif0;
    tif1 = tif1.*Mask;
    intensity(n) = sum(tif1(:));
end
clear tif0 tif1

% 7. plot x
X = [1:1:row]';
figure('color','w');
plot(X, intensity);
xlim([0, 2050]);
xlabel('Frames','fontsize',10);
ylabel('Intensity','fontsize',10);

%%
c = length(find(Mask(:)~=0));
average_intensity = intensity/c;
figure('color','w');
plot(X, average_intensity);
xlim([0, 2050]);
xlabel('Frames','fontsize',10);
ylabel('Intensity','fontsize',10);

%% -- Laplace and iLaplace dROI for Current info

Current = intensity2current(intensity, row);
% clear Intensity imgNum
disp('***''Average each dROI'' has finished***');


%% -- plot Current
prompt = 'Please input the beginning voltage:\n ';
BeginVolt = input(prompt);
prompt = 'Please input the middle voltage:\n ';
EndVolt = input(prompt);
Voltage  = calculateVolt(Current, BeginVolt, EndVolt); % calaulate the X axis - Voltage

figure
plot(Voltage, Current);
title('Graph of current calculated by SPR intensity'); % plot title
xlabel('Voltage/V') % x-axis label
ylabel('Current/A') % y-axis label
disp('***''plot Current'' has finished***');

