tic

Fs = 106;
result  = 'E:\20190306_MoS2_1123_CH18SH\result\';

%%
[~, begin.pike] = max(diff(data(:, 1)));
begin.CS1 = 30770;
begin.start = ceil((begin.CS1 - begin.pike + 1)/10000*106);
begin.end = begin.start +1024 -1;
A3.data = data;
A3.begin = begin;


%% import image Sequence

A3.beginFrame = begin.start;
A3.endFrame = begin.end;
A3.frame = (1 : (A3.endFrame-A3.beginFrame+1))';

[~, A3.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
A3.tifDir = dir(fullfile(A3.tifFile, '*.tiff'));
tif0 = double(imread(fullfile(A3.tifFile, A3.tifDir(1).name)));

%% Total partly sampling analysis vs. bg
A3.roiMask = 'I:\20190306_MoS2_1123_CH18SH\Mask_B5';
[~, A3.roiNames] = ReadTifFileNames(A3.roiMask);

A3.bgMask = 'I:\20190306_MoS2_1123_CH18SH\Mask_B5_bg';
[~, A3.bgNames] = ReadTifFileNames(A3.bgMask);

A3.tifSeq = cell(size(A3.roiNames));
tic
for n = 1:length(A3.roiNames)
    roiMask = ~imread(fullfile(A3.roiMask, A3.roiNames{n}));
    bgMask = ~imread(fullfile(A3.bgMask, A3.bgNames{n}));
    
    check = xor(roiMask, bgMask);
    if sum(check(:)) == 0
        return
    end
    clear check
    
    for ii = A3.beginFrame:A3.endFrame
        tif  = double(imread(fullfile(A3.tifFile, A3.tifDir(ii).name))) - tif0;
        A3.tifSeq{n, 1}((ii-A3.beginFrame+1), 1) = ROImean(tif, roiMask);
        A3.tifSeq{n, 1}((ii-A3.beginFrame+1), 2) = ROImean(tif, bgMask);
    end
    
    [f1, P1] = fft_P1(A3.tifSeq{n, 1}(:, 1), Fs);
    figure('color', 'w');
    plot(f1, P1);
    hold on
    [f2, P2] = fft_P1(A3.tifSeq{n, 1}(:, 2), Fs);
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
    
    A3.sampleArea = cell(size(A3.frame));
    for ii = A3.beginFrame:A3.endFrame
        tif  = double(imread(fullfile(A3.tifFile, A3.tifDir(ii).name))) - tif0;
        A3.sampleArea{(ii-A3.beginFrame) + 1, 1} = ii - A3.beginFrame + 1;
        A3.sampleArea{(ii-A3.beginFrame) + 1, 2} = tif((row(1):row(2)), (col(1):col(2)));
    end
    
end

    
%% Try adpmedian
n = 20;

img = A3.sampleAreaFFT_f{n};
% subplot(1, 2, 1);
subplot(1, 3, 1);
imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
set(gca, 'CLim', [-60 60]);
colormap jet
% map = colormap('jet');
% colorbar;
% imshow(img, 'InitialMagnification', 'fit', 'colormap', hsv);
impixelinfo

img1 = A3.sampleAreaFFT_f{n+1};
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
img = A3.sampleArea{ii, 2};
figure('color', 'w');
imshow(img)
h_img = imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
img1 = fft2(img);
img2 = fftshift(img1);
figure('color', 'w');
imshow(img2, 'DisplayRange', [], 'InitialMagnification', 'fit');
h = drawrectangle;
A3.sampleMask = createMask(h, h_img);


%% Total FFT

tic
sampleArea = A3.sampleArea;
sampleMask = A3.sampleMask;
sampleAreaFFT = cell(size(sampleArea, 1), 1);
parfor n = 1:size(sampleArea, 1)
    sampleAreaFFT{n, 1} = FFTconvert(sampleArea{n, 2}, sampleMask);
end
A3.sampleAreaFFT =  sampleAreaFFT;
clear  sampleAreaFFT  sampleArea
toc

%% Figure

n = 91;
img = A3.sampleArea{n, 2};
subplot(2, 1, 1);
imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
set(gca, 'CLim', [-60 60]);
colormap jet
impixelinfo

m = n - 6;
img1 = A3.sampleAreaFFT_f{m, 1};
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

A3.sampleAreaF1 = cell(size(A3.sampleArea));
for n = 1:size(A3.sampleArea, 1)
    img = A3.sampleArea{n, 2};
    A3.sampleAreaF1{n, 1} = n;
    A3.sampleAreaF1{n, 2} = adpmedian(img, 3);
end

%% Try lowpass Analysis
for n = 1:length(A3.tifSeq)
    temp = A3.tifSeq{n, 1};
    roi = lowp(temp(:, 1), 1, 36, 0.1, 20, Fs);
    bg = lowp(temp(:, 2), 1, 36, 0.1, 20, Fs);
    figure('color', 'w');
    plot(A3.frame, roi, A3.frame, bg)
    legend('roi', 'bg')
    title(n)
    xlabel('Frames')
    ylabel('Intensity')
end

%% Try Time-Frequency Analysis

temp = A3.tifSeq{n, 1};
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

a = A3.sampleAreaF1{1,2};
a(:,:,2) = A3.sampleAreaF1{2,2};

A = {zeros(10),ones(10)};
reshape(cell2mat(A),[size(A{1}) numel(A)]);

%% Total Time-Frequency Analysis - fault

temp = A3.sampleAreaF1(:, 2);
A3.sampleAreaF1_3D = reshape(cell2mat(temp),[size(temp{1}) numel(temp)]);
clear temp

tic
A3.sampleAreaF1_3Df = lowp_s(A3.sampleAreaF1_3D, Fs);
toc

temp = A3.sampleAreaFFT;
A3.sampleAreaF1_3D = reshape(cell2mat(temp),[size(temp{1}) numel(temp)]);
clear temp
A3.sampleAreaFFT1 = lowp_s(A3.sampleAreaF1_3D, Fs);

%% Total Time-Frequency Analysis

% temp = B2.sampleAreaFFT(6:1029);
temp = A3.sampleAreaFFT;
% sampleAreaFFT_3D = zeros([size(B2.sampleMask) 1024]);
for n = 1:length(temp)
    A3.sampleAreaFFT_3D(:, :, n) = temp{n, 1};
end
clear temp

row = size(A3.sampleAreaFFT_3D, 3);
A3.sampleAreaFFT_3D_slide = zeros(row, 2);
for n = 1:row
    A3.sampleAreaFFT_3D_slide(n, 1) = ROImean(A3.sampleAreaFFT_3D(:, :, n), BW1);
    A3.sampleAreaFFT_3D_slide(n, 2) = ROImean(A3.sampleAreaFFT_3D(:, :, n), BW2);
end

tic
A3.sampleAreaFFT_3Df = lowp_s(A3.sampleAreaFFT_3D, Fs);
toc

tic
temp1 = A3.sampleAreaFFT_3Df;
sampleAreaFFT_f = cell(size(A3.sampleAreaFFT));
parfor n = 1:row
    sampleAreaFFT_f{n, 1} = temp1(:, :, n);
end
A3.sampleAreaFFT_f = sampleAreaFFT_f;
clear temp1 sampleAreaFFT_f
toc

%% Total differential intensity analysis 1
BW1 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL1.tif');
BW2 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL2.tif');
BW3 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL3.tif');
row = size(A3.sampleAreaFFT_f, 1);
A3.sampleAreaFFT_fL1L2 = zeros(row, 3);
for n = 1:row
    A3.sampleAreaFFT_fL1L2(n, 1) = ROImean(A3.sampleAreaFFT_f{n}, BW1);
    A3.sampleAreaFFT_fL1L2(n, 2) = ROImean(A3.sampleAreaFFT_f{n}, BW2);
    A3.sampleAreaFFT_fL1L2(n, 3) = ROImean(A3.sampleAreaFFT_f{n}, BW3);
end
L1_L2 = A3.sampleAreaFFT_fL1L2(:, 1) - A3.sampleAreaFFT_fL1L2(:, 2);
L2_L3 = A3.sampleAreaFFT_fL1L2(:, 2) - A3.sampleAreaFFT_fL1L2(:, 3);
L1_L3 = A3.sampleAreaFFT_fL1L2(:, 1) - A3.sampleAreaFFT_fL1L2(:, 3);
x1 = ceil(A3.begin.pike+(A3.beginFrame/106*10000)-1);
x2 = ceil(A3.begin.pike+(A3.beginFrame/106*10000)-1+1023/106*10000);
A3.t1 = (x1/10000:0.0001:x2/10000)';
A3.t1_y = A3.data(x1:x2, 2);
x3 = (A3.beginFrame:1:A3.endFrame)';
A3.t2 = zeros(1024, 1);
for ii = 1:1024
    A3.t2(ii) = (A3.begin.pike+x3(ii)/106*10000-1)/10000;
end
figure('color', 'w');
subplot(2, 1, 1);
plot(A3.t1, 10*A3.t1_y, '-', A3.t2, L1_L2, A3.t2, L2_L3, A3.t2, L1_L3)
subplot(2, 1, 2);
plot(A3.t1, 10*A3.t1_y, '-', A3.t2, A3.sampleAreaFFT_fL1L2(:, 1), A3.t2, A3.sampleAreaFFT_fL1L2(:, 2), A3.t2, A3.sampleAreaFFT_fL1L2(:, 3))

%% Try differential intensity analysis 2
A3.sampleAreaFFT_f_diff = cell(size(A3.sampleAreaFFT_f, 1)-1,1);
for ii = 1:(1024-1)
    A3.sampleAreaFFT_f_diff{ii} = A3.sampleAreaFFT_f{ii+1} - A3.sampleAreaFFT_f{ii};
end

%% Total get Phase difference of Voltage vs time and Intensity vs time
n = 2048;
skip = length(A3.t1)/n;
A3.t1_y0 = A3.t1_y(1:skip:end);
A3.PhDiff = phdiffmeasure(A3.t1_y0(1:1024), L1_L2);

%% Try phase difference 1
[Amp2, Ph2] = freqrespmeasure(B2.t1_y0(1:1024), B2.L1_L2);
[A3.Amp, A3.Ph] = freqrespmeasure(t1_y0(1:1024), A3.L1_L2);


%% Try phase difference 2
tic
output = freqrespmeasure_s(A3.sampleArea_3D, A3.t1_y0);
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
A3.phaseDiff = output;
toc


%%
toc
disp(['Running time: ',num2str(toc)]);

MailToMe('nona1588@outlook.com');


%% Total the experiment of Ru(III)
tifFile = 'E:\20181116_MoS2_CH18-Au\D3_Z1_Ru_PBS_0 -0-4V_0-1VpS_2c_HMMT_200fps';
[A3.tifFolder, A3.tifNames] = ReadTifFileNames(tifFile);
A3.tif0 = double(imread(fullfile(A3.tifFolder, A3.tifNames{1})));
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
m = 2;
A3.beginFrame = 417;
A3.endFrame = 2017;
% C1.sampleArea = cell(C1.endFrame - C1.beginFrame + 1,1);
A3.sampleArea = cell(1600/20+1,2);
[row, col] = ImageJroiLocation(sROI{m});
n = 1;
for ii = A3.beginFrame:20:A3.endFrame
    tif  = double(imread(fullfile(A3.tifFolder, A3.tifNames{ii}))) - A3.tif0;
    A3.sampleArea{n, 1} = n;
    A3.sampleArea{n, 2} = tif((row(1):row(2)), (col(1):col(2)));
    n = n + 1;
end
%%
D3_C1.sampleArea = cell(81, 1);
for ii = 1:81
    D3_C1.sampleArea{ii} = A3.sampleArea{ii, 2} - C1.sampleArea{ii, 2};
end