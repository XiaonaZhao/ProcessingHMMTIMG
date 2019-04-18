tic

Fs = 106;
result  = 'E:\20190306_MoS2_1123_CH18SH\result\';

%%
[~, begin.pike] = max(diff(data(:, 1)));
begin.CS1 = 30770;
begin.start = ceil((begin.CS1 - begin.pike + 1)/10000*106);
begin.end = begin.start +1024 -1;
B5.data = data;
B5.begin = begin;


%% import image Sequence

B5.beginFrame = begin.start;
B5.endFrame = begin.end;
B5.frame = (1 : (B5.endFrame-B5.beginFrame+1))';

[~, B5.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
B5.tifDir = dir(fullfile(B5.tifFile, '*.tiff'));
tif0 = double(imread(fullfile(B5.tifFile, B5.tifDir(1).name)));

%% Total partly sampling analysis vs. bg
B5.roiMask = 'I:\20190306_MoS2_1123_CH18SH\Mask_B5';
[~, B5.roiNames] = ReadTifFileNames(B5.roiMask);

B5.bgMask = 'I:\20190306_MoS2_1123_CH18SH\Mask_B5_bg';
[~, B5.bgNames] = ReadTifFileNames(B5.bgMask);

B5.tifSeq = cell(size(B5.roiNames));
tic
for n = 1:length(B5.roiNames)
    roiMask = ~imread(fullfile(B5.roiMask, B5.roiNames{n}));
    bgMask = ~imread(fullfile(B5.bgMask, B5.bgNames{n}));
    
    check = xor(roiMask, bgMask);
    if sum(check(:)) == 0
        return
    end
    clear check
    
    for ii = B5.beginFrame:B5.endFrame
        tif  = double(imread(fullfile(B5.tifFile, B5.tifDir(ii).name))) - tif0;
        B5.tifSeq{n, 1}((ii-B5.beginFrame+1), 1) = ROImean(tif, roiMask);
        B5.tifSeq{n, 1}((ii-B5.beginFrame+1), 2) = ROImean(tif, bgMask);
    end
    
    [f1, P1] = fft_P1(B5.tifSeq{n, 1}(:, 1), Fs);
    figure('color', 'w');
    plot(f1, P1);
    hold on
    [f2, P2] = fft_P1(B5.tifSeq{n, 1}(:, 2), Fs);
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
    
    B5.sampleArea = cell(size(B5.frame));
    for ii = B5.beginFrame:B5.endFrame
        tif  = double(imread(fullfile(B5.tifFile, B5.tifDir(ii).name))) - tif0;
        B5.sampleArea{(ii-B5.beginFrame) + 1, 1} = ii - B5.beginFrame + 1;
        B5.sampleArea{(ii-B5.beginFrame) + 1, 2} = tif((row(1):row(2)), (col(1):col(2)));
    end
    
end

    
%% Try adpmedian
n = 20;

img = B5.sampleAreaFFT_f{n};
% subplot(1, 2, 1);
subplot(1, 3, 1);
imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
set(gca, 'CLim', [-60 60]);
colormap jet
% map = colormap('jet');
% colorbar;
% imshow(img, 'InitialMagnification', 'fit', 'colormap', hsv);
impixelinfo

img1 = B5.sampleAreaFFT_f{n+1};
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
img = B5.sampleArea{ii, 2};
figure('color', 'w');
imshow(img)
h_img = imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
img1 = fft2(img);
img2 = fftshift(img1);
figure('color', 'w');
imshow(img2, 'DisplayRange', [], 'InitialMagnification', 'fit');
h = drawrectangle;
B5.sampleMask = createMask(h, h_img);


%% Total FFT

tic
sampleArea = B5.sampleArea;
sampleMask = B5.sampleMask;
sampleAreaFFT = cell(size(sampleArea, 1), 1);
parfor n = 1:size(sampleArea, 1)
    sampleAreaFFT{n, 1} = FFTconvert(sampleArea{n, 2}, sampleMask);
end
B5.sampleAreaFFT =  sampleAreaFFT;
clear  sampleAreaFFT  sampleArea
toc

%% Figure

n = 91;
img = B5.sampleArea{n, 2};
subplot(2, 1, 1);
imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
set(gca, 'CLim', [-60 60]);
colormap jet
impixelinfo

m = n - 6;
img1 = B5.sampleAreaFFT_f{m, 1};
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

B5.sampleAreaF1 = cell(size(B5.sampleArea));
for n = 1:size(B5.sampleArea, 1)
    img = B5.sampleArea{n, 2};
    B5.sampleAreaF1{n, 1} = n;
    B5.sampleAreaF1{n, 2} = adpmedian(img, 3);
end

%% Try lowpass Analysis
for n = 1:length(B5.tifSeq)
    temp = B5.tifSeq{n, 1};
    roi = lowp(temp(:, 1), 1, 36, 0.1, 20, Fs);
    bg = lowp(temp(:, 2), 1, 36, 0.1, 20, Fs);
    figure('color', 'w');
    plot(B5.frame, roi, B5.frame, bg)
    legend('roi', 'bg')
    title(n)
    xlabel('Frames')
    ylabel('Intensity')
end

%% Try Time-Frequency Analysis

temp = B5.tifSeq{n, 1};
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

a = B5.sampleAreaF1{1,2};
a(:,:,2) = B5.sampleAreaF1{2,2};

A = {zeros(10),ones(10)};
reshape(cell2mat(A),[size(A{1}) numel(A)]);

%% Total Time-Frequency Analysis - fault

temp = B5.sampleAreaF1(:, 2);
B5.sampleAreaF1_3D = reshape(cell2mat(temp),[size(temp{1}) numel(temp)]);
clear temp

tic
B5.sampleAreaF1_3Df = lowp_s(B5.sampleAreaF1_3D, Fs);
toc

temp = B5.sampleAreaFFT;
B5.sampleAreaF1_3D = reshape(cell2mat(temp),[size(temp{1}) numel(temp)]);
clear temp
B5.sampleAreaFFT1 = lowp_s(B5.sampleAreaF1_3D, Fs);

%% Total Time-Frequency Analysis

% temp = B2.sampleAreaFFT(6:1029);
temp = B5.sampleAreaFFT;
% sampleAreaFFT_3D = zeros([size(B2.sampleMask) 1024]);
for n = 1:length(temp)
    B5.sampleAreaFFT_3D(:, :, n) = temp{n, 1};
end
clear temp

row = size(B5.sampleAreaFFT_3D, 3);
B5.sampleAreaFFT_3D_slide = zeros(row, 2);
for n = 1:row
    B5.sampleAreaFFT_3D_slide(n, 1) = ROImean(B5.sampleAreaFFT_3D(:, :, n), BW1);
    B5.sampleAreaFFT_3D_slide(n, 2) = ROImean(B5.sampleAreaFFT_3D(:, :, n), BW2);
end

tic
B5.sampleAreaFFT_3Df = lowp_s(B5.sampleAreaFFT_3D, Fs);
toc

tic
temp1 = B5.sampleAreaFFT_3Df;
sampleAreaFFT_f = cell(size(B5.sampleAreaFFT));
parfor n = 1:row
    sampleAreaFFT_f{n, 1} = temp1(:, :, n);
end
B5.sampleAreaFFT_f = sampleAreaFFT_f;
clear temp1 sampleAreaFFT_f
toc

%% Total differential intensity analysis 1
BW1 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL1.tif');
BW2 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL2.tif');
BW3 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL3.tif');
row = size(B5.sampleAreaFFT_f, 1);
B5.sampleAreaFFT_fL1L2 = zeros(row, 3);
for n = 1:row
    B5.sampleAreaFFT_fL1L2(n, 1) = ROImean(B5.sampleAreaFFT_f{n}, BW1);
    B5.sampleAreaFFT_fL1L2(n, 2) = ROImean(B5.sampleAreaFFT_f{n}, BW2);
    B5.sampleAreaFFT_fL1L2(n, 3) = ROImean(B5.sampleAreaFFT_f{n}, BW3);
end
L1_L2 = B5.sampleAreaFFT_fL1L2(:, 1) - B5.sampleAreaFFT_fL1L2(:, 2);
L2_L3 = B5.sampleAreaFFT_fL1L2(:, 2) - B5.sampleAreaFFT_fL1L2(:, 3);
L1_L3 = B5.sampleAreaFFT_fL1L2(:, 1) - B5.sampleAreaFFT_fL1L2(:, 3);
x1 = ceil(B5.begin.pike+(B5.beginFrame/106*10000)-1);
x2 = ceil(B5.begin.pike+(B5.beginFrame/106*10000)-1+1023/106*10000);
B5.t1 = (x1/10000:0.0001:x2/10000)';
B5.t1_y = B5.data(x1:x2, 2);
x3 = (B5.beginFrame:1:B5.endFrame)';
B5.t2 = zeros(1024, 1);
for ii = 1:1024
    B5.t2(ii) = (B5.begin.pike+x3(ii)/106*10000-1)/10000;
end
figure('color', 'w');
subplot(2, 1, 1);
plot(B5.t1, 10*B5.t1_y, '-', B5.t2, L1_L2, B5.t2, L2_L3, B5.t2, L1_L3)
subplot(2, 1, 2);
plot(B5.t1, 10*B5.t1_y, '-', B5.t2, B5.sampleAreaFFT_fL1L2(:, 1), B5.t2, B5.sampleAreaFFT_fL1L2(:, 2), B5.t2, B5.sampleAreaFFT_fL1L2(:, 3))

%% Try differential intensity analysis 2
B5.sampleAreaFFT_f_diff = cell(size(B5.sampleAreaFFT_f, 1)-1,1);
for ii = 1:(1024-1)
    B5.sampleAreaFFT_f_diff{ii} = B5.sampleAreaFFT_f{ii+1} - B5.sampleAreaFFT_f{ii};
end

%% Total get Phase difference of Voltage vs time and Intensity vs time
skip = length(B5.t1)/1024;
B5.t1_y0 = B5.t1_y(1:skip:end);
B5.PhDiff = phdiffmeasure(B5.t1_y0(1:1024), L1_L2);

%% Try phase difference 1
[Amp2, Ph2] = freqrespmeasure(B2.t1_y0(1:1024), B2.L1_L2);
[B5.Amp, B5.Ph] = freqrespmeasure(t1_y0(1:1024), B5.L1_L2);


%% Try phase difference 2
tic
output = freqrespmeasure_s(B5.sampleArea_3D, B5.t1_y0);
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
B5.phaseDiff = output;
toc


%%
toc
disp(['Running time: ',num2str(toc)]);

MailToMe('nona1588@outlook.com');


%% Total the experiment of Ru(III)
tifFile = 'E:\20181116_MoS2_CH18-Au\D3_Z1_Ru_PBS_0 -0-4V_0-1VpS_2c_HMMT_200fps';
[B5.tifFolder, B5.tifNames] = ReadTifFileNames(tifFile);
B5.tif0 = double(imread(fullfile(B5.tifFolder, B5.tifNames{1})));
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
m = 2;
B5.beginFrame = 417;
B5.endFrame = 2017;
% C1.sampleArea = cell(C1.endFrame - C1.beginFrame + 1,1);
B5.sampleArea = cell(1600/20+1,2);
[row, col] = ImageJroiLocation(sROI{m});
n = 1;
for ii = B5.beginFrame:20:B5.endFrame
    tif  = double(imread(fullfile(B5.tifFolder, B5.tifNames{ii}))) - B5.tif0;
    B5.sampleArea{n, 1} = n;
    B5.sampleArea{n, 2} = tif((row(1):row(2)), (col(1):col(2)));
    n = n + 1;
end
%%
D3_C1.sampleArea = cell(81, 1);
for ii = 1:81
    D3_C1.sampleArea{ii} = B5.sampleArea{ii, 2} - C1.sampleArea{ii, 2};
end