function Ampl = edl(expName, tifPath, sROI, begin, saveRoute)


L = 1024;

Value.tifFile = tifPath;
Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));
Value.validDir = Value.tifDir(begin.frame:(begin.frame + L -1));


tif0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
for ii = begin.frame:(begin.frame + L -1)
    temp = double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name)));
    Value.tif{(ii-begin.frame+1), 1} = (temp - tif0)./tif0;
end
clear temp

tif_A_peak = zeros(size(tif0, 1), size(tif0, 2));

tif_3D = zeros(size(tif0, 1), size(tif0, 2), L);
for ii = 1:L
    tif_3D(:, :, ii) = Value.tif{ii, 1};
end

for ii = 1:size(tif0, 1)
    parfor jj = 1:size(tif0, 2)
        X = reshape(tif_3D(ii, jj, :), L, 1);
        Y = fft(X);
        P2 = abs(Y/L);
        tif_A_peak(ii, jj) = max(2*P2(2:end-1));
    end
end
Value.A_peak = tif_A_peak;
clear tif_3D

img = figure('color', 'w');
imshow(tif_A_peak, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet
figPath = [saveRoute '\' expName];
saveas(img, figPath, 'fig')
figName = [saveRoute '\' expName '.tif'];
imwrite(tif_A_peak, figName);

Ampl = zeros(length(sROI), 1);
maxAmpl = cell(length(sROI), 1);
for mm = 1:length(sROI)
    [row, col] = ImageJroiLocation(sROI{mm});
    
    tifValue = zeros(L, 1);
    for ii = begin.frame:(begin.frame + L -1)
        temp = Value.tif{(ii-begin.frame+1), 1};
        temp = temp((row(1):row(2)), (col(1):col(2)));
        tifValue(ii, 1) = mean(temp(:));
    end
    clear temp temp1
    Y = fft(tifValue);
    P2 = abs(Y/L);
    Ampl(mm, 1) = max(2*P2(2:end-1));
    
    temp = tif_A_peak((row(1):row(2)), (col(1):col(2)));
    maxAmpl{mm, 1} = temp(:);
    
end
Value.ampl = Ampl;
Value.maxAmpl = maxAmpl;

cellpath = [saveRoute '\' expName '.mat']; 
save(cellpath, 'Value', '-v7.3');

close all


end



% 
% % function MoS2()
% 
% % I want to have fft filter for every frame
% tic
% 
% hwait = waitbar(0, 'Testing!!!!!!!! >>>>>>>>');
% 
% Fs = 100;
% 
% % get TIF
% [~, B5.tifFile] = uigetfile('*.tif', '*.tiff', 'Multiselect', 'on', 'Read tif Folder');
% B5.tifDir = dir(fullfile(B5.tifFile, '*.tif'));
% B5.tif0 = double(imread(fullfile(B5.tifFile, B5.tifDir(1).name)));
% B5.beginFrame = 674;
% B5.endFrame = 2274;
% B5.frame = (1 : (B5.endFrame-B5.beginFrame+1))';
% 
% [~, B6.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
% B6.tifDir = dir(fullfile(B6.tifFile, '*.tif'));
% B6.tif0 = double(imread(fullfile(B6.tifFile, B6.tifDir(1).name)));
% B6.beginFrame = 174;
% B6.endFrame = 1774;
% B6.frame = (1 : (B6.endFrame-B6.beginFrame+1))';
% 
% % [cstrFilenames, cstrPathname] = uigetfile(...
% %     {'*.*',  'All Files (*.*)';...
% %     '*.zip',  'Zip-files (*.zip)';...
% %     '*.roi',  'ROI (*.roi)'...
% %     },'Pick a .roi imageJ file');
% % [sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
% % 
% % [B5.row, B5.col] = ImageJroiLocation(sROI{1});
% % [B6.row, B6.col] = ImageJroiLocation(sROI{2});
% 
% % mask = ~imread('F:\MoS2_final\MoS2_20190513_CH18S-Au\Mask\MaskB6.tif');
% % Bgroup.mask = mask;
% 
% sample = cell(1601, 1);
% for ii = 1:1601
%     temp = B5.beginFrame + ii - 1;
%     tif1 = double(imread(fullfile(B5.tifFile, B5.tifDir(temp).name))) - B5.tif0;
% %     sampleArea1 = tif((B5.row(1):B5.row(2)), (B5.col(1):B5.col(2)));
%     
%     temp = B6.beginFrame + ii - 1;
%     tif2 = double(imread(fullfile(B6.tifFile, B6.tifDir(temp).name))) - B6.tif0;
% %     sampleArea2 = tif((B6.row(1):B6.row(2)), (B6.col(1):B6.col(2)));
%     
%     sample{ii, 1} = tif2 - tif1;
% end
% 
% % % get FFT MASK
% % 
% % 
% % % get FFT filter
% % sample  = Bgroup.sample;
% % sampleFFT = cell(size(sample));
% % parfor ii = 1:size(sample, 1)
% %     sampleFFT{ii, 1} = FFTconvert(sample{ii, 1}, fftMask);
% % end
% % clear sample
% 
% sampleFFT = sample;
% clear sample
% 
% %  get array reshape
% sampleFFT_3D = zeros(1024, 1024, 1601);
% parfor ii = 1:length(temp)
%     sampleFFT_3D(:, :, ii) = imboxfilt(sampleFFT{ii, 1}, 11);
% end
% clear sampleFFT
% 
% % low pass filter in timeline
% sampleFFT_3Df = lowp_s(sampleFFT_3D, Fs);
% clear sampleFFT_3D
% 
% 
% % reshape timeline to spaceslide
% sampleFFT_f = cell(size(sampleFFT_3Df, 1), 1);
% parfor ii = 1:size(sampleFFT_3Df, 1)
%     sampleFFT_f{ii, 1} = sampleFFT_3Df(:, :, ii);
% end
% clear sampleFFT_3Df
% 
% 
% 
% save('sampleFFT_f.mat', 'sampleFFT_f', '-v7.3');
% 
% close(hwait)
% toc