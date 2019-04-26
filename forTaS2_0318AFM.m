prompt = {'Reference number of Data in Exp group:',...
    'Enter the folder of Masks for the Data:',...
    'Enter the matlab file of Data timeline:',...
    'Enter the scan rate of cyclic voltammetry (mV):',...
    'Save route:'};

dlg_title = 'Input';
num_lines = 1;

defaultans = {'A1',...
    'F:\TaS2\20190324_TaS2_0318_ITO\MaskA1',...
    'G:\TaS2\20190318_TaS2_ITO\matlab_data\B2.mat',...
    '100',...
    'F:\TaS2\20190324_TaS2_0318_ITO\result'};

answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
expName = answer{1};
maskPath = answer{2};
loader = answer{3};
rate = str2num(answer{4});
saveRoute = answer{5};

% if DC
TaS2_DC(expName, maskPath, saveRoute);
% if CV
% TaS2(expName, maskPath, loader, rate, saveRoute);

%%
load('F:\TaS2\20190324_TaS2_0318_ITO\result\expTab_TaS2_0324.mat')

for m = 1:size(expTab_TaS2, 1)
    expName = expTab_TaS2{m, 2};
    tifPath = expTab_TaS2{m, 3};
    maskPath = expTab_TaS2{m, 4};
    saveRoute = expTab_TaS2{m, 5};
    TaS2_DC(expName, tifPath, maskPath, saveRoute);
end

%%
[~, A3.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
A3.tifDir = dir(fullfile(A3.tifFile, '*.tiff'));
tif0 = double(imread(fullfile(A3.tifFile, A3.tifDir(1).name)));

[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
m =1;
[row, col] = ImageJroiLocation(sROI{m});

A3.beginFrame1 = 1001;
A3.endFrame1 = 1350;
A3.part1 = cell(350, 1);
for ii = A3.beginFrame1:A3.endFrame1
    tif  = double(imread(fullfile(A3.tifFile, A3.tifDir(ii).name))) - tif0;
    A3.part1{(ii-A3.beginFrame1) + 1, 1} = ii - A3.beginFrame1 + 1;
    A3.part1{(ii-A3.beginFrame1) + 1, 2} = tif((row(1):row(2)), (col(1):col(2)));
end

A3.part1diff = cell(size(A3.part1, 1) -1, 1);
for ii = 1:length(A3.part1diff)
    A3.part1diff{ii} = A3.part1{ii+1, 2} - A3.part1{ii, 2};
end

temp = A3.part1diff;
A3.part1_3D = zeros([size(img) 349]);
for n = 1:length(temp)
    A3.part1_3D(:, :, n) = temp{n, 1};
end
clear temp

[A3.part1diff_value, A3.part1diff_time] = imgSeqDiff_s(A3.part1_3D);

%%

