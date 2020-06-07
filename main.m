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
saveRoute = 'G:\MoS2\MoS2_0802\_Result_f';
%%
exp = cell(3, 5);
fields = {'expName', 'tifPath', 'tifName', 'sROI', 'begin'};
Value = cell2struct(exp, fields, 2);

%%
Value(1).expName = 'A3';
Value(2).expName = 'A4';
Value(3).expName = 'A5';

%%
% 1. Read .tiff files names
tifFile = 'G:\MoS2\MoS2_0802\A3_Ru_-0-1 -0-3V_100mV_2c_HMMT100fps';
[Value(1).tifPath, Value(1).tifName] = ReadTifFileNames(tifFile);

tifFile = 'G:\MoS2\MoS2_0802\A4_Ru_-0-1 -0-3V_100mV_2c_HMMT100fps';
[Value(2).tifPath, Value(2).tifName] = ReadTifFileNames(tifFile);

tifFile = 'G:\MoS2\MoS2_0802\A5_PBS_-0-1 -0-3V_100mV_2c_HMMT100fps';
[Value(3).tifPath, Value(3).tifName] = ReadTifFileNames(tifFile);

%%
Value(1).begin = 282;
Value(2).begin = 294;
Value(3).begin = 272;

% maskFile = 'E:\20181123_MoS2_CH18\matlabMask';
% [maskFolder, maskNames] = ReadTifFileNames(maskFile);

% if size(tifNames, 1) == size(maskNames, 1)
%     row = size(tifNames, 1);
% else
%     disp('***''The amount of TIFF doesn''t match that of MASK.***');
% end
%%
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));

Value(1).sROI = sROI{1, 4};
Value(2).sROI = sROI{1, 2};
Value(3).sROI = sROI{1, 1};

Value(1).prep = Value(1).begin;
Value(2).prep = Value(3).begin;

%%
tic
L = 850;

[row1, col1] = ImageJroiLocation(Value(1).sROI);
tif = double(imread(fullfile(Value(1).tifPath, Value(1).tifName{1})));
tif01 = tif(row1(1):row1(2), col1(1):col1(2));

[row2, col2] = ImageJroiLocation(Value(2).sROI);
tif = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{1})));
tif02 = tif(row2(1):row2(2), col2(1):col2(2));

tif = double(imread(fullfile(Value(1).tifPath, Value(1).tifName{Value(1).prep})));
base01 = -ROImean(tif(row1(1):row1(2), col1(1):col1(2)), fake);
tif = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{Value(2).prep})));
base02 = -ROImean(tif(row2(1):row2(2), col2(1):col2(2)), fake);

% tif1_0 = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{Value(1).prep})));
% tif2_0 = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{Value(2).prep})));
% tif1_0 = tif1_0((row1(1):row1(2)), (col1(1):col1(2)));
% tif2_0 = tif2_0((row2(1):row2(2)), (col2(1):col2(2)));
% 
% tif1_400 = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{399+Value(1).begin})));
% tif2_400 = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{399+Value(2).begin})));
% tif1_400 = tif1_400((row1(1):row1(2)), (col1(1):col1(2)));
% tif2_400 = tif2_400((row2(1):row2(2)), (col2(1):col2(2)));
% 
% tif1_600 = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{599+Value(1).begin})));
% tif2_600 = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{599+Value(2).begin})));
% tif1_600 = tif1_600((row1(1):row1(2)), (col1(1):col1(2)));
% tif2_600 = tif2_600((row2(1):row2(2)), (col2(1):col2(2)));

exp = cell(L, 1);
% for ii = 1:L
%     tif = double(imread(fullfile(Value(1).tifPath, Value(1).tifName{ii+Value(1).begin-1})));
%     tif1 = tif(row1(1):row1(2), col1(1):col1(2)) - tif01;
%     base1 = -ROImean(tif1, fake);
%     tif = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{ii+Value(2).begin-1})));
%     tif2 = tif(row2(1):row2(2), col2(1):col2(2))- tif02;
%     base2 = -ROImean(tif2, fake);
%     
%     exp{ii, 1} = tif1 - (tif2*(base01/base02).*ifake0 + tif2*(base1/base2).*fake0);
% end

for ii = 1:L
    tif = double(imread(fullfile(Value(1).tifPath, Value(1).tifName{ii+Value(1).begin-1})));
    tif1 = tif(row1(1):row1(2), col1(1):col1(2)) - tif01;
    
    tif = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{ii+Value(2).begin-1})));
    tif2 = tif(row2(1):row2(2), col2(1):col2(2))- tif02;
    
    exp{ii, 1} = tif1 - 3*tif2;
end


cellpath = [saveRoute '\' Value(1).expName '_' Value(2).expName '.mat'];
save(cellpath, 'exp', '-v7.3');
toc

%%
X = (1:L)';
Qc = zeros(L, 1);
parfor ii = 1:L
    temp = exp{ii, 1};
    Qc(ii, 1) = sum(temp(:));
end
Ia = diff(Qc);

figure('color','white');
plot(X, Qc)
xlim([0 850])
xlabel('Frames')
ylabel('Qc (a.u.)')

figure('color','white');
plot(potential, Qc(1:801, 1))
xlabel('Potential (V)')
ylabel('Qc (a.u.)')

figure('color','white');
plot(X(2:end), Ia)
xlim([0 850])
xlabel('Frames')
ylabel('Ia (a.u.)')
%%
% figure('color','white');
Ib = lowp(Ia, 16*0.1, 15, 0.01, 20, 100);
plot(potential(401:800), Ib(401:800))
xlabel('Potential (V)')
ylabel('Ia (a.u.)')


%%
tic
exp_3D = zeros(size(tif1, 1), size(tif1, 2), L);
parfor ii = 1:L
    exp_3D(:, :, ii) = exp{ii, 1};
end
toc

tic
Current_3D = zeros(size(tif1, 1), size(tif1, 2), L-1);
Current_3D_lowp = zeros(size(tif1, 1), size(tif1, 2), L-1);
for ii = 1:size(tif1, 1)
    parfor jj = 1:size(tif1, 2)
        Intensity = reshape(exp_3D(ii, jj, :), L, 1);
        Current_3D(ii, jj, :) = intensity2current(Intensity, L);
        Current_3D_lowp(ii, jj, :) = lowp(Current_3D(ii, jj, :), 1*0.1, 6, 0.01, 10, 100);
    end
end
clear exp_3D
toc

tic
Current = cell(L-1, 1);
Current_lowp = cell(L-1, 1);
parfor ii = 1:(L-1)
    Current{ii, 1} = -reshape(Current_3D(:, :, ii), [size(tif1, 1), size(tif1, 2)]);
    Current_lowp{ii, 1} = -reshape(Current_3D_lowp(:, :, ii), [size(tif1, 1), size(tif1, 2)]);
end
toc

tic
cellpath = [saveRoute '\' Value(1).expName '_' Value(2).expName '_Current.mat'];
save(cellpath, 'Current', '-v7.3');
cellpath = [saveRoute '\' Value(1).expName '_' Value(2).expName '_Current_lowp.mat'];
save(cellpath, 'Current_lowp', '-v7.3');
clear Current_3D Current_3D_lowp
toc


%%
tic
%%
TifSaveRoute = 'G:\MoS2\MoS2_0802\_Result_f\TIF_A3';
parfor ii = 1:length(exp)
    c = num2str(ii, '%06d');  
    TifName = [TifSaveRoute '\Inte_', c, '.tif'];
    I = exp{ii, 1};
%     I = I/255;
    temp = I - ones(size(I))*min(I(:));
    I = uint8(temp/max(temp(:))*255);
    imwrite(I, TifName);
end
%%
TifSaveRoute = 'G:\MoS2\MoS2_0802\_Result\TIF_B3_B4_Current';
parfor ii = 1:length(Current)
    c = num2str(ii, '%06d');
    TifName = [TifSaveRoute '\Cur_', c, '.tif'];
    A = Current{ii, 1};
    A = A/255;
    imwrite(A, TifName);
end
toc

%%
ii = 600;
I = exp{ii, 1};
A = Current{ii, 1};
B = Current_lowp{ii, 1};
% A(A>200) = 200;
% A(A<0) = 0;
B(B>1000) = 1000;
B(B<-1500) = -1500;
figure('color','white');
% subplot(131)
imshow(A, 'DisplayRange',[], 'InitialMagnification', 'fit');
colormap jet
colorbar
% lim_A = caxis;
caxis([-2500 1500])
% caxis([-1500 1000])

%%

figure('color','white');
subplot(151)
imshow(Current{436, 1}, 'DisplayRange',[], 'InitialMagnification', 'fit');
colormap gray
caxis([-2500 1500])

subplot(152)
imshow(Current{559, 1}, 'DisplayRange',[], 'InitialMagnification', 'fit');
colormap gray
caxis([-2500 1500])

subplot(153)
imshow(Current{635, 1}, 'DisplayRange',[], 'InitialMagnification', 'fit');
colormap gray
caxis([-2500 1500])

subplot(154)
imshow(Current{765, 1}, 'DisplayRange',[], 'InitialMagnification', 'fit');
colormap gray
caxis([-2500 1500])

subplot(155)
imshow(Current{836, 1}, 'DisplayRange',[], 'InitialMagnification', 'fit');
colormap gray
caxis([-2500 1500])

%%
% figure('color','white');


    figure('color','white');
for ii = 1:5

    subplot(2,5,ii)
    temp = img_plot{ii, 1};
%     temp = temp.*temp.*temp;
    imshow(temp, 'DisplayRange',[], 'InitialMagnification', 'fit');
    colormap parula(5)
    caxis([-1500 1000])
end

    figure('color','white');
for ii = 6:10

    subplot(2,5,ii)
    temp = img_plot{ii, 1};
%     temp = temp.*temp;
    imshow(temp, 'DisplayRange',[], 'InitialMagnification', 'fit');
    colormap parula(5)
    caxis([-2500 1500])
end

%%
for ii = 11:15
    figure('color','white');
%     subplot(3,5,ii)
    temp = img_plot{ii, 1};
    temp = temp.*temp.*temp;
    imshow(temp, 'DisplayRange',[], 'InitialMagnification', 'fit');
    colormap default
    caxis([-900^2.5 900^2.5])
end



%%
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));

%%
Intensity_line = zeros(length(exp), length(sROI));
for ii = 1:length(exp)
    temp  = exp{ii, 1};
    for jj = 1:length(sROI)
        [row, col] = ImageJroiLocation(sROI{1, jj});
        ROI = temp((row(1):row(2)), (col(1):col(2)));
        Intensity_line(ii, jj) = mean(ROI(:));
    end
end
%%
Current_line = zeros(length(exp)-1, length(sROI));
for ii = 1:length(sROI)
    Current_line(:, ii) = -intensity2current(Intensity_line(:, ii), L);
end

%%
d = 19;
% figure('color','white');
plot(potential(401:end, 1), 0.01*Current_line(400:800, length(sROI)), '.',...
    potential(401:d:end, 1), 0.01*Current_line(401:d:800, 1), 'o',...
    CHI_V, 45000*CHI_A);
xlabel('Potential (V)')
ylabel('Ia (a.u.)')
% current_lowp_line = zeros(length(Current), length(sROI));
% for ii = 1:length(Current)
%     for jj = 1:length(sROI)
%         [row, col] = ImageJroiLocation(sROI{1, jj});
%         temp = Current{ii, 1}((row(1):row(2)), (col(1):col(2)));
%         current_line(ii, jj) = mean(temp(:));
%         temp = Current_lowp{ii, 1}((row(1):row(2)), (col(1):col(2)));
%         current_lowp_line(ii, jj) = mean(temp(:));
%     end
% end

%%
rate = 100;
Fs = 100;
potential = potentialLine(rate, Fs, -0.1, -0.3); % 2c
%%
X = (1:length(Current))';
figure('color','white');
% d = 17;
% plot(potential(401:d:end, 1), current_lowp_line(436:d:836, 4), 'o');
hold on
for ii = 1:length(sROI)
% ii = 1;
%     line = 0.001*lowp(current_line(:, ii), 16*0.1, 15, 0.1, 20, 100);
    line = 0.001*lowp(current_line(:, ii), 1*0.1, 6, 0.01, 10, 100);
    plot(potential(401:end, 1), line(400:800, 1));
%     plot(potential(401:end, 1), current_line(400:800, ii));
%     plot(X(400:end, 1), current_line(400:end, ii));
%     plot(X(400:end, 1), current_lowp_line(400:end, ii));
end
xlabel('Potential (V)')
ylabel('I2Current (a.u.)')
hold off


%%
figure('color','white');
hold on
for ii = 1:length(sROI)
    plot(potential(399:end, 1), current_lowp_line(399:end, ii));
end
hold off
%%
deBaseline = zeros(size(Intensity_line));
for ii = [1 3 4 10 12 13]
    line = lowp(Intensity_line(:, ii), 40*0.1, 15, 0.1, 20, 100);
    deBaseline(:, ii) = driftBaseline(X, line);
end
%%
figure('color','white');
hold off
% for ii = 1:length(sROI)
for ii = [1 3 4 10 12 13]
%     line = lowp(Intensity_line(:, ii), 40*0.1, 15, 0.1, 20, 100);
    plot(X, deBaseline(:, ii));
    hold on
end

%%
figure('color', 'white');
hold on
% for ii = 1:length(sROI)
    ii = 4;
    line = lowp(current_line(:, ii), 1*0.1, 6, 0.01, 10, 100);
%     line = smooth(line, 38);
    plot(X(399:end, 1), current_line(399:end, ii))
    plot(X(399:end, 1), line(399:end, 1))
%     plot(potential(400:800, 1), line(400:800, 1))
    % d = 22;
    % plot(X(399:d:end, 1), line(399:d:end, 1), 'o')
% end
hold off
%%
for ii = 1:length(sROI)
    [f1, P1] = fft_P1(current_line(:, ii), 100);
    figure('color', 'white');
    plot(f1, P1);
    xlim([0 50])
    title(['Spectrum of ROI ' num2str(ii)])
    xlabel('f (Hz)','fontsize',10)
    ylabel('|P(f)|','fontsize',10)
end
%%
mask = ~imread('G:\MoS2\MoS2_0802\_ROI\Mask_f.tif');
Curve = zeros(L, 1);
parfor ii = 1:L
    Curve(ii, 1) = -ROImean(exp{ii, 1}, mask);
end

%%
% fX = (1:L)';
% figure('color', 'white');
A4f = lowp(A3, 15*0.1, 15, 0.01, 20, 100);
A4f = smooth(A4f, 31);
% plot(X, Curve, X, filterCurve)
plot(X, A3, X, A4f)
%%
ii = 523;
tif_ii = exp{ii, 1} - exp{399, 1};
redCon_ii = (-tif_ii)*10/(-(filterCurve(600)-filterCurve(400)));
localSums = imboxfilt(redCon_ii, 17);
localSums0 = localSums;
localSums = abs(localSums);
localSums(localSums > 10) = 10;
imshow(localSums, 'DisplayRange',[], 'InitialMagnification', 'fit');
colormap jet

%%
savepath = 'G:\MoS2\MoS2_0802\_Result\TIF_B3_B4_concentration_529x509\';
for ii = 802:1602
%     ii = 600;
%     tif_ii = exp{ii, 1};
    tif_ii = exp{ii, 1} - exp{399, 1};
    % tif1_n = B1_A4.sampleArea_cut{n};
    % tif1_n = double(imread(fullfile(tifFolder, tifNames{n}))) - tif0;
    % % tif_n = double(imread(fullfile(tifFolder, tifNames{n-2})));
    %
    % mask = double(imread(fullfile(maskFolder, maskNames{1}))/255);
    % tif1_n = tif1_n.*mask;
    
    if ii > 600
        Voltage = -0.3 + ((ii-600)/100)*0.1;
    else
        Voltage = -0.1 - ((ii-400)/100)*0.1;
    end
    
%     tif_ii(tif_ii > 0) = 0;
%     m = -500;
%     tif_ii(tif_ii < m) = m;
    
    redCon_ii = (-tif_ii/2)*10/(-(filterCurve(600)-filterCurve(400)));
    localSums = imboxfilt(redCon_ii, 11);
%     localSums0 = localSums;
%     localSums = abs(localSums);
    localSums(localSums > 10) = 10;
%     figure('color', 'w');
    imshow(localSums, 'DisplayRange',[], 'InitialMagnification', 'fit');
    title([num2str(Voltage), ' V']);
    colormap jet
    map=colormap('jet');
    colorbar;
%     impixelinfo
    set(gca, 'CLim', [-300 -40]);
    h=colorbar;
    set(get(h,'title'),'string','Ru(II), mM', 'FontSize', 12);
%     pause(0.05);
%     saveas(gcf,[savepath, 'B3_B4_' num2str(ii, '%04d'), '.tif']);
end

%%
video = VideoWriter('demo.avi'); %初始化一个avi文件
video.FrameRate = 24;
open(video);
for ii = 1:401  %图像序列个数
    frame = imread(fullfile(tifFolder, tifNames{ii}));
    writeVideo(video, frame);
end
close(video);


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
    tif01 = double(imread(fullfile(fileFolder, fileNames{n}))) - tif0;
    tif01 = tif01.*Mask;
    intensity(n) = sum(tif01(:));
end
clear tif0 tif1

% 7. plot x
X = (1:row)';
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
Current_3D = zeros(size(Y(2:end,:)));
for n = 1:col
    % Current(:,n) = intensity2current(intensity(:,n), row);
    Current_3D(:,n) = intensity2current(filterY(:,n), 801);
    % clear Intensity imgNum
end
disp('***''Average each dROI'' has finished***');



%% -- plot Current
prompt = 'Please input the beginning voltage:\n ';
BeginVolt = input(prompt);
prompt = 'Please input the middle voltage:\n ';
EndVolt = input(prompt);
Voltage  = calculateVolt(Current_3D, BeginVolt, EndVolt); % calaulate the X axis - Voltage

%%
figure('color','w');
% plot(Voltage, Current);
for n = 1:col
    %     plot(X, Y(:, n)); % get raw lines
    plot(X(4:end), -Current_3D((3:end), n), '.'); % get fitted lines
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
