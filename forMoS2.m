tic

Fs = 106;
result  = 'H:\20190306_MoS2_1123_CH18SH\result\';

%% import image Sequence
B2.beginFrame = 160;
B2.endFrame = 1400;
B2.frame = (1 : (B2.endFrame-B2.beginFrame+1))';

[~, B2.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
B2.tifDir = dir(fullfile(B2.tifFile, '*.tiff'));
tif0 = double(imread(fullfile(B2.tifFile, B2.tifDir(1).name)));

%% Total partly sampling analysis vs. bg
B2.roiMask = 'H:\20190306_MoS2_1123_CH18SH\Mask_B2';
[~, B2.roiNames] = ReadTifFileNames(B2.roiMask);

B2.bgMask = 'H:\20190306_MoS2_1123_CH18SH\Mask_B2_bg';
[~, B2.bgNames] = ReadTifFileNames(B2.bgMask);

B2.tifSeq = cell(size(B2.roiNames));
for n = 1:length(B2.roiNames)
    roiMask = ~imread(fullfile(B2.roiMask, B2.roiNames{n}));
    bgMask = ~imread(fullfile(B2.bgMask, B2.bgNames{n}));
    
    check = xor(roiMask, bgMask);
    if sum(check(:)) == 0
        return
    end
    clear check
    
    for ii = B2.beginFrame:B2.endFrame
        tif  = double(imread(fullfile(B2.tifFile, B2.tifDir(ii).name))) - tif0;
        B2.tifSeq{n, 1}((ii-B2.beginFrame+1), 1) = ROImean(tif, roiMask);
        B2.tifSeq{n, 1}((ii-B2.beginFrame+1), 2) = ROImean(tif, bgMask);
    end
    
    [f1, P1] = fft_P1(B2.tifSeq{n, 1}(:, 1), Fs);
    figure('color', 'w');
    plot(f1, P1);
    hold on
    [f2, P2] = fft_P1(B2.tifSeq{n, 1}(:, 2), Fs);
    plot(f2, P2);
    xlim([0, 20]); ylim([0, 50]);
    xlabel('f (Hz)','fontsize',10)
    ylabel('|P(f)|','fontsize',10)
    set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2)
    set(gca, 'linewidth', 1.5)
    title(n)
    legend('ROI', 'Background')
    hold off
    
end


%% Total sampleArea
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
for m = 1:length(sROI)
    [row, col] = ImageJroiLocation(sROI{m});
    
    B2.sampleArea = cell(size(B2.frame) - 1);
    for ii = (B2.beginFrame+1):B2.endFrame
        tif  = double(imread(fullfile(B2.tifFile, B2.tifDir(ii).name))) - tif0;
        B2.sampleArea{(ii-B2.beginFrame), 1} = ii - B2.beginFrame;
        B2.sampleArea{(ii-B2.beginFrame), 2} = tif((row(1):row(2)), (col(1):col(2)));
    end
    
end

%% Try adpmedian

img = B2.sampleArea{95, 2};
subplot(1, 2, 1);
imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet
% map = colormap('jet');
% colorbar;
% imshow(img, 'InitialMagnification', 'fit', 'colormap', hsv);
impixelinfo

img1 = adpmedian(img, 3);
subplot(1, 2, 2);
imshow(img1, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet
impixelinfo
%% Total FFT
sampleArea = B2.sampleArea;
sampleMask = B2.sampleMask;
sampleAreaFFT = cell(size(sampleArea, 1), 1);
parfor n = 1:size(sampleArea, 1)
    sampleAreaFFT{n, 1} = FFTconvert(sampleArea{n, 2}, sampleMask);
end
B2. sampleAreaFFT =  sampleAreaFFT;
clear  sampleAreaFFT  sampleArea

%%
n = 316;
img = B2.sampleArea{n, 2};
subplot(1, 2, 1);
imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet
impixelinfo

img1 = B2.sampleAreaFFT_f{n, 1};
% img1 = B2.sampleAreaF1_3Df(:,:,n);
% img1 = B2.sampleAreaF1_3Df2(:,:,n);
subplot(1, 2, 2);
imshow(img1,'DisplayRange', [],'InitialMagnification', 'fit')
% imshow(img1, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet
impixelinfo

%% Total adpmedian

B2.sampleAreaF1 = cell(size(B2.sampleArea));
for n = 1:size(B2.sampleArea, 1)
    img = B2.sampleArea{n, 2};
    B2.sampleAreaF1{n, 1} = n;
    B2.sampleAreaF1{n, 2} = adpmedian(img, 3);
end

%% Try lowpass Analysis
for n = 1:length(B2.tifSeq)
    temp = B2.tifSeq{n, 1};
    roi = lowp(temp(:, 1), 1, 36, 0.1, 20, Fs);
    bg = lowp(temp(:, 2), 1, 36, 0.1, 20, Fs);
    figure('color', 'w');
    plot(B2.frame, roi, B2.frame, bg)
    legend('roi', 'bg')
    title(n)
    xlabel('Frames')
    ylabel('Intensity')
end

%% Try Time-Frequency Analysis

% temp = B2.tifSeq{n, 1};
% roi = lowp(temp(:, 1), 4, 12, 0.1, 20, Fs);
roi = lowp(temp(:, 1), 4, 16, 33*0.01, 20, Fs);
roi = highpass(roi, 5, Fs);
bg = lowp(temp(:, 2), 4, 16, 33*0.01, 20, Fs);
% bg = lowp(temp(:, 2), 3, 7, 0.1, 20, Fs);
bg = highpass(bg, 5, Fs);
% figure('color', 'w');
plot(x, roi, x, bg)
legend('roi', 'bg')
title(n)
xlim([0 1000])
xlabel('Frames')
ylabel('Intensity')
% figure('color', 'w');spectrogram(roi,[],[],[],Fs);
% figure('color', 'w');spectrogram(bg,[],[],[],Fs);

%% Try one-dot Time-Frequency Analysis
% from 2D to 3D

a = B2.sampleAreaF1{1,2};
a(:,:,2) = B2.sampleAreaF1{2,2};

A = {zeros(10),ones(10)};
reshape(cell2mat(A),[size(A{1}) numel(A)]);

%% Total Time-Frequency Analysis - fault

temp = B2.sampleAreaF1(:, 2);
B2.sampleAreaF1_3D = reshape(cell2mat(temp),[size(temp{1}) numel(temp)]);
clear temp

B2.sampleAreaF1_3Df = lowp_s(B2.sampleAreaF1_3D, Fs);

temp = B2.sampleAreaFFT;
B2.sampleAreaF1_3D = reshape(cell2mat(temp),[size(temp{1}) numel(temp)]);
clear temp
B2.sampleAreaFFT1 = lowp_s(B2.sampleAreaF1_3D, Fs);

%% Total Time-Frequency Analysis

temp = B2.sampleAreaFFT(6:1029);
% sampleAreaFFT_3D = zeros([size(B2.sampleMask) 1024]);
for n = 1:length(temp)
    B2.sampleAreaFFT_3D(:, :, n) = temp{n, 1};
end
clear temp

row = size(B2.sampleAreaFFT_3D, 3);
B2.sampleAreaFFT_3D_slide = zeros(row, 2);
for n = 1:row
    B2.sampleAreaFFT_3D_slide(n, 1) = ROImean(B2.sampleAreaFFT_3D(:, :, n), BW1);
     B2.sampleAreaFFT_3D_slide(n, 2) = ROImean(B2.sampleAreaFFT_3D(:, :, n), BW2);
end

B2.sampleAreaFFT_3Df = lowp_s(B2.sampleAreaFFT_3D, Fs);

temp1 = B2.sampleAreaFFT_3Df;
for n = 1:row
    B2.sampleAreaFFT_f{n, 1} = temp1(:, :, n);
end
clear temp1

%%
toc
disp(['Running time: ',num2str(toc)]);

MailToMe('nona1588@outlook.com');