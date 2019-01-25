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
tifFile = 'E:\20181227_MoS2_CH18-SH\recount_IMG\mono_-0-4V';
tifFolder = fullfile(tifFile);
dirOutput = dir(fullfile(tifFolder, '*.tif'));
dirFileTif = sortObj(dirOutput);
tifNames = {dirFileTif.name}';

maskFile = 'E:\20181227_MoS2_CH18-SH\recount_IMG\mask-mono';
maskFolder = fullfile(maskFile);
dirOutput = dir(fullfile(maskFolder, '*.tif'));
dirFileMask = sortObj(dirOutput);
maskNames = {dirFileMask.name}';

if size(tifNames, 1) == size(maskNames, 1)
    row = size(tifNames, 1);
else
    disp('***''The amount of TIFF doesn''t match that of MASK.***');
end

%%
% 2. Read one .tiff file
% 3. Use a fixed mask on the image into new one
% 4. Sum the the new image to return x
% 5. Release the original image
% 6. loop 2 to 5
intensity = zeros(row, 1);
Mask =~ imread('F:\20190116_MoS2_CH18-SH_0108\MaskA1.tif');
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
% xlim([0, 2050]);
xlabel('Frames','fontsize',10);
ylabel('Intensity','fontsize',10);

%% Get average
c = length(find(Mask(:)~=0));
average_intensity = intensity/c;
figure('color','w');
plot(X, average_intensity);
xlim([0, 2050]);
xlabel('Frames','fontsize',10);
ylabel('Intensity','fontsize',10);

%% substract sampling CV with empty CV, performed only at Wenli Lv's PC
% Input sampling CV
fileFolder1 = fullfile('F:\20190116_MoS2_CH18-SH_0108\A1_z1_10mMRu250mMPBNa_0 -0-4V_0-1VpS_2c_MoS2_CH18-S-Au_HMMT1600fps');
dirOutput1 = dir(fullfile(fileFolder1,'*.tif'));
dirFile1 = sortObj(dirOutput1);
fileNames1 = {dirFile1.name}';
tifBegin1 = double(imread(fullfile(fileFolder1, fileNames1{1})));
% Input empty background CV
fileFolder2 = fullfile('F:\20190116_MoS2_CH18-SH_0108\B1_empty_10mMRu250mMPBNa_0 -0-4V_0-1VpS_2c_MoS2_CH18-S-Au_HMMT1600fps');
dirOutput2 = dir(fullfile(fileFolder2,'*.tif'));
dirFile2 = sortObj(dirOutput2);
fileNames2 = {dirFile2.name}';
tifBegin2 = double(imread(fullfile(fileFolder2, fileNames2{1})));

row = 16*1600 + 1;
path = 'F:\20190116_MoS2_CH18-SH_0108\tiff_A1subB1';
for n = 1:row
    tifSample1 = double(imread(fullfile(fileFolder1, fileNames1{n+4258}))) - tifBegin1; % uniformization
    tifSample2 = double(imread(fullfile(fileFolder2, fileNames2{n+2945}))) - tifBegin2; % uniformization
    tifSample = tifSample1 - tifSample2; % matlab cannot process 32-bit TIFF files.
    file = fileNames1{n};
    pathfile = fullfile(path, file);
    imwrite(tifSample, pathfile, 'tiff');
end



%% -- Laplace and iLaplace dROI for Current info
Current = zeros(size(Y(2:end,:)));
for n = 1:col
    % Current(:,n) = intensity2current(intensity(:,n), row);
    Current(:,n) = intensity2current(filterY(:,n), 801);
    % clear Intensity imgNum
end
disp('***''Average each dROI'' has finished***');



%% -- plot Current
prompt = 'Please input the beginning voltage:\n ';
BeginVolt = input(prompt);
prompt = 'Please input the middle voltage:\n ';
EndVolt = input(prompt);
Voltage  = calculateVolt(Current, BeginVolt, EndVolt); % calaulate the X axis - Voltage
%%
figure('color','w');
% plot(Voltage, Current);
for n = 1:col
    %     plot(X, Y(:, n)); % get raw lines
    plot(X(4:end), -Current((3:end), n), '.'); % get fitted lines
    hold on
end
title('Graph of current calculated by SPR intensity'); % plot title
xlabel('Voltage/V') % x-axis label
ylabel('Current/A') % y-axis label
disp('***''plot Current'' has finished***');

%%
s = zeros(row, 1);
for n = 1:row
    img = imread(fullfile(tifFolder, tifNames{n}));
    mask = imread(fullfile(maskFolder, maskNames{n}));
    s(n) = calculateROIintensity(img, mask);
end
