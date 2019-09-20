% function MoS2()

% I want to have fft filter for every frame
tic

hwait = waitbar(0, 'Testing!!!!!!!! >>>>>>>>');

Fs = 100;

% get TIF
[~, B5.tifFile] = uigetfile('*.tif', '*.tiff', 'Multiselect', 'on', 'Read tif Folder');
B5.tifDir = dir(fullfile(B5.tifFile, '*.tif'));
B5.tif0 = double(imread(fullfile(B5.tifFile, B5.tifDir(1).name)));
B5.beginFrame = 674;
B5.endFrame = 2274;
B5.frame = (1 : (B5.endFrame-B5.beginFrame+1))';

[~, B6.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
B6.tifDir = dir(fullfile(B6.tifFile, '*.tif'));
B6.tif0 = double(imread(fullfile(B6.tifFile, B6.tifDir(1).name)));
B6.beginFrame = 174;
B6.endFrame = 1774;
B6.frame = (1 : (B6.endFrame-B6.beginFrame+1))';

% [cstrFilenames, cstrPathname] = uigetfile(...
%     {'*.*',  'All Files (*.*)';...
%     '*.zip',  'Zip-files (*.zip)';...
%     '*.roi',  'ROI (*.roi)'...
%     },'Pick a .roi imageJ file');
% [sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
% 
% [B5.row, B5.col] = ImageJroiLocation(sROI{1});
% [B6.row, B6.col] = ImageJroiLocation(sROI{2});

% mask = ~imread('F:\MoS2_final\MoS2_20190513_CH18S-Au\Mask\MaskB6.tif');
% Bgroup.mask = mask;

sample = cell(1601, 1);
for ii = 1:1601
    temp = B5.beginFrame + ii - 1;
    tif1 = double(imread(fullfile(B5.tifFile, B5.tifDir(temp).name))) - B5.tif0;
%     sampleArea1 = tif((B5.row(1):B5.row(2)), (B5.col(1):B5.col(2)));
    
    temp = B6.beginFrame + ii - 1;
    tif2 = double(imread(fullfile(B6.tifFile, B6.tifDir(temp).name))) - B6.tif0;
%     sampleArea2 = tif((B6.row(1):B6.row(2)), (B6.col(1):B6.col(2)));
    
    sample{ii, 1} = tif2 - tif1;
end

% % get FFT MASK
% 
% 
% % get FFT filter
% sample  = Bgroup.sample;
% sampleFFT = cell(size(sample));
% parfor ii = 1:size(sample, 1)
%     sampleFFT{ii, 1} = FFTconvert(sample{ii, 1}, fftMask);
% end
% clear sample

sampleFFT = sample;
clear sample

%  get array reshape
sampleFFT_3D = zeros(1024, 1024, 1601);
parfor ii = 1:length(temp)
    sampleFFT_3D(:, :, ii) = imboxfilt(sampleFFT{ii, 1}, 11);
end
clear sampleFFT

% low pass filter in timeline
sampleFFT_3Df = lowp_s(sampleFFT_3D, Fs);
clear sampleFFT_3D


% reshape timeline to spaceslide
sampleFFT_f = cell(size(sampleFFT_3Df, 1), 1);
parfor ii = 1:size(sampleFFT_3Df, 1)
    sampleFFT_f{ii, 1} = sampleFFT_3Df(:, :, ii);
end
clear sampleFFT_3Df



save('sampleFFT_f.mat', 'sampleFFT_f', '-v7.3');

close(hwait)
toc