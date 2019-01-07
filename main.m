% main function for Intensity2Current
clc
clear

% ATTENTION:

% There is still two problems in this programme. one is that the CURRENT
% calculated by the iLaplace doesn't match with the result from the
% electrochemical workstation. The other is that the present background
% cutting doesn't correct the bias from the charge and discharge of gold
% plate.

% Author: nonazhao@mail.ustc.edu.cn;
% Created: 29 June 2018

% -- Important modification: 22 Nov 2018 --
% added submain.m, fft_P1.m, lowp.m and several plot lines


%% -- Read the alive .tiff files

prompt = 'Increment of frames for the samples in alive period:\n ';
skipNum = input(prompt);
[imgSeq, imgNum] = ReadTifFiles(...
    'Open sampling image sequence', skipNum); % uint16 cell


% -- Remove the no-electro background

[BgSeq, BgNum] = ReadTifFiles('Open background sequence', 0); % uint16 cell
imgSubtractBg = BgdRemoval(imgSeq, imgNum, BgSeq, BgNum);
clear imgSeq BgSeq BgNum skipNum
disp('***''Read the alive .tiff files'' has finished***');


%% -- Select ROIs

[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.roi',  'ROI (*.roi)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.*',  'All Files (*.*)',...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
switch sROI.strType
    case 'Rectangle'
        RectBounds = sROI.vnRectBounds;
        % [x-left_up y-left_up x-right_down y-right_down]
        col = [RectBounds(2), RectBounds(2), RectBounds(4), RectBounds(4)];
        row = [RectBounds(1), RectBounds(3), RectBounds(3), RectBounds(1)];
    case 'Polygon'
        Polygon = sROI.mnCoordinates;
        col = Polygon(:,2);
        row = Polygon(:,1);
end
BW = roipoly(imgSubtractBg{1}, col, row);
nonzeroBW  = length(find(BW(:)~=0));
BW = BW*1;
% "*1" turns logical into double, then "uint16" turn double into uint16.
% set more variable to monitor the matrix changes
imgSegment = cell(imgNum,1);
for j = 1:imgNum
    imgSegment{j} = imgSubtractBg{j}.*BW; % double cell
end
clear cstrFilenames cstrPathname sROI
clear col row BW j imgSubtractBg
disp('***''Select ROIs'' has finished***');


%% -- Average each dROI

Intensity = averROI(imgSegment, imgNum, nonzeroBW);
clear imgSegment nonzeroBW
disp('***''Select ROIs'' has finished***');


%% -- Laplace and iLaplace dROI for Current info

Current = intensity2current(Intensity, imgNum);
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

