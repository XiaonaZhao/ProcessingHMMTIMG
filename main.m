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
fileFolder = fullfile('D:\Guaduate\SPR-others\Experiments\20190108\B3_N3_10mMRu250mMPBNa_0 -0-4V_2c_0-1VpS_MoS2_CH18-S-Au_sp80_HMMT100fps');
dirOutput = dir(fullfile(fileFolder,'*.tif'));
fileNames = {dirOutput.name}';
row = size(fileNames, 1);

% 2. Read one .tiff file
% 3. Use a fixed mask on the image into new one
% 4. Sum the the new image to return x
% 5. Release the original image
% 6. loop 2 to 5
intensity = zeros(row, 1);
Mask =~ imread('D:\Guaduate\SPR-others\Experiments\20190108\mask_0108_B3.tif');
for n = 1: row
    tif0 = double(imread(fullfile(fileFolder, fileNames{n})));
    tif1 = tif0.*Mask;
    intensity(n) = sum(tif1(:));
end
clear tif0 tif1

% 7. plot x
X = [1:1:row]';
figure('color','w');
plot(X, intensity);


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

