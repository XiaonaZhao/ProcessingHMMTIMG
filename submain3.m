%% This page works for the Figure 2 and the SI video
% clear all
load('G:\MoS2\MoS2_0802\_Result_f\matlab_drifted_fake.mat') % monolayer
%%
% B3_Ru_-0-1 -0-3V_100mV_2c_HMMT100fps
[row1, col1] = ImageJroiLocation(Value(1).sROI);
tif = double(imread(fullfile(Value(1).tifPath, Value(1).tifName{1})));
tif01 = tif(row1(1):row1(2), col1(1):col1(2));

% B4_PBS_-0-1 -0-3V_100mV_2c_HMMT100fps
[row2, col2] = ImageJroiLocation(Value(2).sROI);
tif = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{1})));
tif02 = tif(row2(1):row2(2), col2(1):col2(2));

tif = double(imread(fullfile(Value(1).tifPath, Value(1).tifName{Value(1).prep})));
base01 = -ROImean(tif(row1(1):row1(2), col1(1):col1(2)), fake0);
tif = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{Value(2).prep})));
base02 = -ROImean(tif(row2(1):row2(2), col2(1):col2(2)), fake0);

exp = cell(L, 1);
for ii = 1:L
    tif = double(imread(fullfile(Value(1).tifPath, Value(1).tifName{ii+Value(1).begin-1})));
    tif1 = tif(row1(1):row1(2), col1(1):col1(2)) - tif01;
    
    tif = double(imread(fullfile(Value(2).tifPath, Value(2).tifName{ii+Value(2).begin-1})));
    tif2 = tif(row2(1):row2(2), col2(1):col2(2))- tif02;
    
%     exp{ii, 1}= (tif1 - 3*tif2)./tif01;
    exp{ii, 1}= tif1 - 3.5*tif2;
end
%% Inport a smaller roi
% G:\MoS2\MoS2_0802\_ROI\RoiSet.zip
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
[row, col] = ImageJroiLocation(sROI{1, 6}); % Or roi8 monolayers
% [row, col] = ImageJroiLocation(sROI{1, 4}); % A3 multilayers


%% preparations for imaging the intensity
figure('color', 'w');
ii = 600;
tif_ii = exp{ii, 1};
tif_ii = tif_ii(row(1):row(2), col(1):col(2));
if ii > 600
    Voltage = -0.3 + ((ii-600)/100)*0.1;
else
    Voltage = -0.1 - ((ii-400)/100)*0.1;
end
localSums = imboxfilt(tif_ii, 13);
localSums(localSums > 0) = 0;
imshow(localSums, 'DisplayRange',[], 'InitialMagnification', 'fit');
title([num2str(Voltage), ' V']);
c = fire; % turbo parula fire
c = flipud(c);
map = colormap(c);
colorbar;
%     impixelinfo
% set(gca, 'CLim', [-1200 0]);
h = colorbar;
set(get(h,'title'),'string','Ru(II)/mM', 'FontSize', 12);


%% preparations for imaging the concentration
ii = 540;
% tif_ii = exp{ii, 1};
tif_ii = exp{ii, 1} - exp{801, 1}; % It is quite important
tif_ii = tif_ii(row(1):row(2), col(1):col(2))/2;

if ii > 600
    Voltage = -0.3 + ((ii-600)/100)*0.1;
else
    Voltage = -0.1 - ((ii-400)/100)*0.1;
end

% redCon_ii = (-tif_ii)*10/600;
redCon_ii = (-tif_ii)*10/(-(filterCurve(600)-filterCurve(400)));
localSums = imboxfilt(redCon_ii, 13);
localSums(localSums > 10) = 10;
localSums(localSums < 0) = 0;
imshow(localSums, 'DisplayRange',[], 'InitialMagnification', 'fit');
title([num2str(Voltage), ' V'], 'FontSize', 14, 'FontWeight', 'bold');
colormap fire % parula fire
colorbar
set(gca, 'CLim', [0 10]);
h=colorbar;
set(h, 'FontSize', 14, 'FontWeight', 'bold');
set(get(h,'title'),'string','[Ru(NH_3)_6]Cl_3^2^+ (mM)', 'FontSize', 12);


%% Cutscale Vision
savepath = 'G:\MoS2\MoS2_0802\_Result_f\TIF_B3_B4_fire_RuII\';
for ii = 400:800
% ii = 600;
    tif_ii = exp{ii, 1} - exp{801, 1};
    tif_ii = tif_ii(row(1):row(2), col(1):col(2))/2;
    
    if ii > 600
        Voltage = -0.3 + ((ii-600)/100)*0.1;
    else
        Voltage = -0.1 - ((ii-400)/100)*0.1;
    end
    
    redCon_ii = (-tif_ii)*10/(-(filterCurve(600)-filterCurve(400)));
    localSums = imboxfilt(redCon_ii, 13);
    
    localSums(localSums > 10) = 10;
    localSums(localSums < 0) = 0;
    
    imshow(localSums, 'DisplayRange',[], 'InitialMagnification', 'fit');
    title(['t = ' num2str((ii-400)/100), ' s, Voltage = ' num2str(Voltage), ' V'], 'FontSize', 14, 'FontWeight', 'bold');
    colormap fire % parula
    colorbar;

    set(gca, 'CLim', [0 10]);
    h = colorbar;
    set(h, 'FontSize', 14, 'FontWeight', 'bold');
    set(get(h,'title'),'string','[Ru(NH_3)_6]Cl_3^2^+ (mM)', 'FontSize', 12);

    pause(0.01);
    saveas(gcf,[savepath, 'B3_B4_' num2str(ii, '%04d'), '.tif']);
    
end

%% Fullscale Vision
savepath = 'G:\MoS2\MoS2_0802\_Result_f\TIF_B3_B4_parula_RuII_fullscale\';
for ii = 400:800

    tif_ii = (exp{ii, 1} - exp{801, 1})/2;
%     tif_ii = tif_ii(row(1):row(2), col(1):col(2))/2;
    
    if ii > 600
        Voltage = -0.3 + ((ii-600)/100)*0.1;
    else
        Voltage = -0.1 - ((ii-400)/100)*0.1;
    end
    
    redCon_ii = (-tif_ii)*10/(-(filterCurve(600)-filterCurve(400)));
    localSums = imboxfilt(redCon_ii, 13);
    
    localSums(localSums > 10) = 10;
    localSums(localSums < 0) = 0;
    
    imshow(localSums, 'DisplayRange',[], 'InitialMagnification', 'fit');
    title(['t = ' num2str((ii-400)/100), ' s, Voltage = ' num2str(Voltage), ' V'], 'FontSize', 14, 'FontWeight', 'bold');
    colormap parula
    colorbar;

    set(gca, 'CLim', [0 10]);
    h = colorbar;
    set(h, 'FontSize', 14, 'FontWeight', 'bold');
    set(get(h,'title'),'string','[Ru(NH_3)_6]Cl_3^2^+ (mM)', 'FontSize', 12);
    scalebar();
    
    pause(0.05);
    saveas(gcf,[savepath, 'B3_B4_' num2str(ii, '%04d'), '.tif']);
    
end

%% Make a video
% tifFolder = 'G:\MoS2\MoS2_0802\_Result_f\TIF_B3_B4_parula_RuII';
tifFolder = 'G:\MoS2\MoS2_0802\_Result_f\TIF_B3_B4_parula_RuII_fullscale';
[~, tifNames] = ReadTifFileNames(tifFolder);
video = VideoWriter('TIF_B3_B4_parula_RuII_fullscale.avi'); % Prepare a .avi file
video.FrameRate = 100; % The same as the fps
open(video);
for ii = 1:401  % The frame number for video
    frame = imread(fullfile(tifFolder, tifNames{ii}));
    writeVideo(video, frame);
end
close(video);
%%
figure('color', 'w');
plot()