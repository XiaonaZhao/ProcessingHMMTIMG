tic

Fs = 106;
result  = 'E:\20190306_MoS2_1123_CH18SH\result\';

%%
[~, begin.pike] = max(diff(data(:, 1)));
begin.CS1 = 30770;
begin.start = ceil((begin.CS1 - begin.pike + 1)/10000*106);
begin.end = begin.start +1024 -1;
C2.data = data;
C2.begin = begin;


%% import image Sequence

C2.beginFrame = begin.start;
C2.endFrame = begin.end;
C2.frame = (1 : (C2.endFrame-C2.beginFrame+1))';

[~, C2.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
C2.tifDir = dir(fullfile(C2.tifFile, '*.tiff'));
tif0 = double(imread(fullfile(C2.tifFile, C2.tifDir(1).name)));

%% Total partly sampling analysis vs. bg
C2.roiMask = 'I:\20190306_MoS2_1123_CH18SH\Mask_B5';
[~, C2.roiNames] = ReadTifFileNames(C2.roiMask);

C2.bgMask = 'I:\20190306_MoS2_1123_CH18SH\Mask_B5_bg';
[~, C2.bgNames] = ReadTifFileNames(C2.bgMask);

C2.tifSeq = cell(size(C2.roiNames));
tic
for n = 1:length(C2.roiNames)
    roiMask = ~imread(fullfile(C2.roiMask, C2.roiNames{n}));
    bgMask = ~imread(fullfile(C2.bgMask, C2.bgNames{n}));
    
    check = xor(roiMask, bgMask);
    if sum(check(:)) == 0
        return
    end
    clear check
    
    for ii = C2.beginFrame:C2.endFrame
        tif  = double(imread(fullfile(C2.tifFile, C2.tifDir(ii).name))) - tif0;
        C2.tifSeq{n, 1}((ii-C2.beginFrame+1), 1) = ROImean(tif, roiMask);
        C2.tifSeq{n, 1}((ii-C2.beginFrame+1), 2) = ROImean(tif, bgMask);
    end
    
    [f1, P1] = fft_P1(C2.tifSeq{n, 1}(:, 1), Fs);
    figure('color', 'w');
    plot(f1, P1);
    hold on
    [f2, P2] = fft_P1(C2.tifSeq{n, 1}(:, 2), Fs);
    plot(f2, P2);
    %     xlim([0, 20]); ylim([0, 50]);% set Fs = 5 Hz
    xlim([0, 30]); ylim([0, 60]);
    xlabel('f (Hz)','fontsize',10)
    ylabel('|P(f)|','fontsize',10)
    set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2)
    set(gca, 'linewidth', 1.5)
    title(n)
    legend('ROI', 'Background')
    hold off
    
end
toc

%% Total sampleArea
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
for m = 1:length(sROI)
    [row, col] = ImageJroiLocation(sROI{m});
    
    C2.sampleArea = cell(size(C2.frame));
    for ii = C2.beginFrame:C2.endFrame
        tif  = double(imread(fullfile(C2.tifFile, C2.tifDir(ii).name))) - tif0;
        C2.sampleArea{(ii-C2.beginFrame) + 1, 1} = ii - C2.beginFrame + 1;
        C2.sampleArea{(ii-C2.beginFrame) + 1, 2} = tif((row(1):row(2)), (col(1):col(2)));
    end
    
end

    
%% Try adpmedian
n = 20;

img = C2.sampleAreaFFT_f{n};
% subplot(1, 2, 1);
subplot(1, 3, 1);
imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
set(gca, 'CLim', [-60 60]);
colormap jet
% map = colormap('jet');
% colorbar;
% imshow(img, 'InitialMagnification', 'fit', 'colormap', hsv);
impixelinfo

img1 = C2.sampleAreaFFT_f{n+1};
% subplot(1, 2, 2);
subplot(1, 3, 2);
imshow(img1, 'DisplayRange', [], 'InitialMagnification', 'fit');
set(gca, 'CLim', [-60 60]);
colormap jet
impixelinfo

subplot(1, 3, 3);
imshow((img1 - img), 'DisplayRange', [], 'InitialMagnification', 'fit');
set(gca, 'CLim', [-15 15]);
colormap jet
impixelinfo

%% Total sampleMask
ii = 44;
img = C2.sampleArea{ii, 2};
figure('color', 'w');
imshow(img)
h_img = imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
img1 = fft2(img);
img2 = fftshift(img1);
figure('color', 'w');
imshow(img2, 'DisplayRange', [], 'InitialMagnification', 'fit');
h = drawrectangle;
C2.sampleMask = createMask(h, h_img);


%% Total FFT

tic
sampleArea = C2.sampleArea;
sampleMask = C2.sampleMask;
sampleAreaFFT = cell(size(sampleArea, 1), 1);
parfor n = 1:size(sampleArea, 1)
    sampleAreaFFT{n, 1} = FFTconvert(sampleArea{n, 2}, sampleMask);
end
C2.sampleAreaFFT =  sampleAreaFFT;
clear  sampleAreaFFT  sampleArea
toc

%% Figure

n = 91;
img = C2.sampleArea{n, 2};
subplot(2, 1, 1);
imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
set(gca, 'CLim', [-60 60]);
colormap jet
impixelinfo

m = n - 6;
img1 = C2.sampleAreaFFT_f{m, 1};
% img1 = B2.sampleAreaF1_3Df(:,:,n);
% img1 = B2.sampleAreaF1_3Df2(:,:,n);
subplot(2, 1, 2);
imshow(img1,'DisplayRange', [],'InitialMagnification', 'fit')
% imshow(img1, 'DisplayRange', [], 'InitialMagnification', 'fit');
set(gca, 'CLim', [-40 40]);
colormap jet
% colorbar
impixelinfo

%% Total adpmedian

C2.sampleAreaF1 = cell(size(C2.sampleArea));
for n = 1:size(C2.sampleArea, 1)
    img = C2.sampleArea{n, 2};
    C2.sampleAreaF1{n, 1} = n;
    C2.sampleAreaF1{n, 2} = adpmedian(img, 3);
end

%% Try lowpass Analysis
for n = 1:length(C2.tifSeq)
    temp = C2.tifSeq{n, 1};
    roi = lowp(temp(:, 1), 1, 36, 0.1, 20, Fs);
    bg = lowp(temp(:, 2), 1, 36, 0.1, 20, Fs);
    figure('color', 'w');
    plot(C2.frame, roi, C2.frame, bg)
    legend('roi', 'bg')
    title(n)
    xlabel('Frames')
    ylabel('Intensity')
end

%% Try Time-Frequency Analysis

temp = C2.tifSeq{n, 1};
% roi = lowp(temp(:, 1), 4, 12, 0.1, 20, Fs);
roi = lowp(temp(:, 1), 4, 16, 33*0.01, 20, Fs);
roi = highpass(roi, 10, Fs);
bg = lowp(temp(:, 2), 4, 16, 33*0.01, 20, Fs);
% bg = lowp(temp(:, 2), 3, 7, 0.1, 20, Fs);
bg = highpass(bg, 10, Fs);
% figure('color', 'w');
plot(x, roi, x, bg)
% plot(x, temp(:, 1), x, temp(:, 2))
legend('roi', 'bg')
title(n)
xlim([200 300])
xlabel('Frames')
ylabel('Intensity')
figure('color', 'w');
subplot(2,1,1)
% spectrogram(temp(:, 1),[],[],[],Fs); title('roi');
spectrogram(roi,[],[],[],Fs); title('roi');
set(gca, 'CLim', [-40 40]);
% colormap jet
% figure('color', 'w');
subplot(2,1,2)
% spectrogram(temp(:, 2),[],[],[],Fs); title('roi');
spectrogram(bg,[],[],[],Fs); title('bg');
set(gca, 'CLim', [-40 40]);
% colormap jet

%% Try one-dot Time-Frequency Analysis
% from 2D to 3D

a = C2.sampleAreaF1{1,2};
a(:,:,2) = C2.sampleAreaF1{2,2};

A = {zeros(10),ones(10)};
reshape(cell2mat(A),[size(A{1}) numel(A)]);

%% Total Time-Frequency Analysis - fault

temp = C2.sampleAreaF1(:, 2);
C2.sampleAreaF1_3D = reshape(cell2mat(temp),[size(temp{1}) numel(temp)]);
clear temp

tic
C2.sampleAreaF1_3Df = lowp_s(C2.sampleAreaF1_3D, Fs);
toc

temp = C2.sampleAreaFFT;
C2.sampleAreaF1_3D = reshape(cell2mat(temp),[size(temp{1}) numel(temp)]);
clear temp
C2.sampleAreaFFT1 = lowp_s(C2.sampleAreaF1_3D, Fs);

%% Total Time-Frequency Analysis

% temp = B2.sampleAreaFFT(6:1029);
temp = C2.sampleAreaFFT;
% sampleAreaFFT_3D = zeros([size(B2.sampleMask) 1024]);
for n = 1:length(temp)
    C2.sampleAreaFFT_3D(:, :, n) = temp{n, 1};
end
clear temp

row = size(C2.sampleAreaFFT_3D, 3);
C2.sampleAreaFFT_3D_slide = zeros(row, 2);
for n = 1:row
    C2.sampleAreaFFT_3D_slide(n, 1) = ROImean(C2.sampleAreaFFT_3D(:, :, n), BW1);
    C2.sampleAreaFFT_3D_slide(n, 2) = ROImean(C2.sampleAreaFFT_3D(:, :, n), BW2);
end

tic
C2.sampleAreaFFT_3Df = lowp_s(C2.sampleAreaFFT_3D, Fs);
toc

tic
temp1 = C2.sampleAreaFFT_3Df;
sampleAreaFFT_f = cell(size(C2.sampleAreaFFT));
parfor n = 1:row
    sampleAreaFFT_f{n, 1} = temp1(:, :, n);
end
C2.sampleAreaFFT_f = sampleAreaFFT_f;
clear temp1 sampleAreaFFT_f
toc

%% Total differential intensity analysis 1
BW1 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL1.tif');
BW2 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL2.tif');
BW3 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL3.tif');
row = size(C2.sampleAreaFFT_f, 1);
C2.sampleAreaFFT_fL1L2 = zeros(row, 3);
for n = 1:row
    C2.sampleAreaFFT_fL1L2(n, 1) = ROImean(C2.sampleAreaFFT_f{n}, BW1);
    C2.sampleAreaFFT_fL1L2(n, 2) = ROImean(C2.sampleAreaFFT_f{n}, BW2);
    C2.sampleAreaFFT_fL1L2(n, 3) = ROImean(C2.sampleAreaFFT_f{n}, BW3);
end
L1_L2 = C2.sampleAreaFFT_fL1L2(:, 1) - C2.sampleAreaFFT_fL1L2(:, 2);
L2_L3 = C2.sampleAreaFFT_fL1L2(:, 2) - C2.sampleAreaFFT_fL1L2(:, 3);
L1_L3 = C2.sampleAreaFFT_fL1L2(:, 1) - C2.sampleAreaFFT_fL1L2(:, 3);
x1 = ceil(C2.begin.pike+(C2.beginFrame/106*10000)-1);
x2 = ceil(C2.begin.pike+(C2.beginFrame/106*10000)-1+1023/106*10000);
C2.t1 = (x1/10000:0.0001:x2/10000)';
C2.t1_y = C2.data(x1:x2, 2);
x3 = (C2.beginFrame:1:C2.endFrame)';
C2.t2 = zeros(1024, 1);
for ii = 1:1024
    C2.t2(ii) = (C2.begin.pike+x3(ii)/106*10000-1)/10000;
end
figure('color', 'w');
subplot(2, 1, 1);
plot(C2.t1, 10*C2.t1_y, '-', C2.t2, L1_L2, C2.t2, L2_L3, C2.t2, L1_L3)
subplot(2, 1, 2);
plot(C2.t1, 10*C2.t1_y, '-', C2.t2, C2.sampleAreaFFT_fL1L2(:, 1), C2.t2, C2.sampleAreaFFT_fL1L2(:, 2), C2.t2, C2.sampleAreaFFT_fL1L2(:, 3))

%% Try differential intensity analysis 2
C2.sampleAreaFFT_f_diff = cell(size(C2.sampleAreaFFT_f, 1)-1,1);
for ii = 1:(1024-1)
    C2.sampleAreaFFT_f_diff{ii} = C2.sampleAreaFFT_f{ii+1} - C2.sampleAreaFFT_f{ii};
end

%% Total get Phase difference of Voltage vs time and Intensity vs time
n = 2048;
skip = length(C2.t1)/n;
C2.t1_y0 = C2.t1_y(1:skip:end);
C2.PhDiff = phdiffmeasure(C2.t1_y0(1:1024), L1_L2);

%% Try phase difference 1
[Amp2, Ph2] = freqrespmeasure(B2.t1_y0(1:1024), B2.L1_L2);
[C2.Amp, C2.Ph] = freqrespmeasure(t1_y0(1:1024), C2.L1_L2);


%% Try phase difference 2
tic
output = freqrespmeasure_s(C2.sampleArea_3D, C2.t1_y0);
% Ma = max(output(:));
% Mi = min(output(:));
% OneMi = ones(size(B5.sampleMask))*Mi;
output = output + (output<0)*2*pi - (output>2*pi)*2*pi;
figure('color', 'w');
imshow(output, 'DisplayRange', [], 'InitialMagnification', 'fit');
% set(gca, 'CLim', [Mi Ma]);
colormap jet
colorbar
impixelinfo
C2.phaseDiff = output;
toc

figure('color','w');
subplot(2, 1, 1); hist(Hz5)
subplot(2, 1, 2); hist(Hz10)

%%
toc
disp(['Running time: ',num2str(toc)]);

MailToMe('nona1588@outlook.com');


%% Total the experiment of Ru(III)
tifFile = 'E:\20181116_MoS2_CH18-Au\D3_Z1_Ru_PBS_0 -0-4V_0-1VpS_2c_HMMT_200fps';
[C2.tifFolder, C2.tifNames] = ReadTifFileNames(tifFile);
C2.tif0 = double(imread(fullfile(C2.tifFolder, C2.tifNames{1})));
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
m = 2;
C2.beginFrame = 417;
C2.endFrame = 2017;
% C1.sampleArea = cell(C1.endFrame - C1.beginFrame + 1,1);
C2.sampleArea = cell(1600/20+1,2);
[row, col] = ImageJroiLocation(sROI{m});
n = 1;
for ii = C2.beginFrame:20:C2.endFrame
    tif  = double(imread(fullfile(C2.tifFolder, C2.tifNames{ii}))) - C2.tif0;
    C2.sampleArea{n, 1} = n;
    C2.sampleArea{n, 2} = tif((row(1):row(2)), (col(1):col(2)));
    n = n + 1;
end
%%
D3_C1.sampleArea = cell(81, 1);
for ii = 1:81
    D3_C1.sampleArea{ii} = C2.sampleArea{ii, 2} - C1.sampleArea{ii, 2};
end