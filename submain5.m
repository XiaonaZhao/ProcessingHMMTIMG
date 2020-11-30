% For Figure
load('G:\MoS2\20181116_MoS2_CH18-Au\TIFF_CV\Cut_D3-C1_mat.mat')
exp = cell(81, 1);
for ii = 1:81
    exp{ii, 1} = D3_C1.sampleArea{ii};
end
%%
load('E:\MoS2\MoS2_0802\_Result_f\matlab_drifted_fake.mat') % monolayer
%%
% figure('color', 'w');
ii = 60;
tif_ii = exp{ii, 1};
if ii > 40
    Voltage = -0.4 + ((ii-41)/10)*0.1;
else
    Voltage = -0.4 - ((ii-41)/10)*0.1;
end
localSums = imboxfilt(tif_ii, 11);
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
%%
ii =  40;
tif_ii = exp{ii, 1};
% tif_ii = exp{ii, 1} - exp{1, 1}; % It is an optional.

if ii > 40
    Voltage = -0.4 + ((ii-41)/10)*0.1;
else
    Voltage = -0.4 - ((ii-41)/10)*0.1;
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
h = colorbar;
set(h, 'FontSize', 14, 'FontWeight', 'bold');
set(get(h,'title'),'string','[Ru(NH_3)_6]Cl_3^2^+ (mM)', 'FontSize', 12);
%%
savepath = 'G:\MoS2\20181116_MoS2_CH18-Au\TIFF_CV\Cut_D3_C1_fire\';
for ii = 1:81
    tif_ii = exp{ii, 1};
    % tif_ii = exp{ii, 1} - exp{1, 1}; % It is an optional.
    
    if ii > 40
        Voltage = -0.4 + ((ii-41)/10)*0.1;
    else
        Voltage = -0.4 - ((ii-41)/10)*0.1;
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
    h = colorbar;
    set(h, 'FontSize', 14, 'FontWeight', 'bold');
    set(get(h,'title'),'string','[Ru(NH_3)_6]Cl_3^2^+ (mM)', 'FontSize', 12);
    
    pause(0.05);
    saveas(gcf,[savepath, 'D3_C1_' num2str(ii, '%04d'), '.tif']);
end
%% D3-C1_f1: ROI1 ROI3 ROI6, intensity
% load('G:\MoS2\20181116_MoS2_CH18-Au\TIFF_CV\D3_C1_f1.mat')
% figure('color', 'w');
Y = lowp(ROI/7000, 3*0.1, 15, 0.01, 20, 100);
plot(X, Y)
xlim([0 1600])
xlabel('Potential (V vs. Ag/AgCl)')
ylabel('\DeltaI/I (a.u.)')
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 1);
set(gca, 'linewidth', 1, 'FontSize', 16)

%% D3-C1_f1: ROI1 ROI3 ROI6, CV
% figure('color', 'w');
Current = zeros(1600, 3);
for ii = 1:3
    Current(:, ii) = -intensity2current(Y(:, ii), 1601);
end
CurrDens = Current/((0.25^2) * (10^(-8))); % Current Density (A/cm2)
% plot(Potential(1:801), CurrDens(800:end,:))
% ylabel('Current density (a.u.)')
plot(Potential(1:801), Current(800:end,:))
xlabel('Potential (V vs. Ag/AgCl)')
ylabel('Current (a.u.)')
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 1);
set(gca, 'linewidth', 1, 'FontSize', 16)
%% D3-C1_f1: ROI6 ROI1 CHI, CV
load('G:\MoS2\20181116_MoS2_CH18-Au\TIFF_CV\D3_C1_f1.mat')
figure('color', 'w');

colororder({...
    '#D95319' ... % [0.8500 0.3250 0.0980]	'#D95319', orange
%     ,'#EDB120'... % '#EDB120' [0.9290 0.6940 0.1250], yellow
    '#000000'... % 'black'	'k'	[0 0 0]	'#000000', black
    % ,'#0072BD'... % '#0072BD' [0 0.4470 0.7410], blue
    }) 
box on

yyaxis left
% plot(Potential(1:801), Current(800:end, 3)) % ROI6, background
hold on
plot(Potential(1:801), Current(800:end, 1)) % ROI1, sample edge
xlim([-0.4 0])
xlabel('Potential (V vs. Ag/AgCl)')
ylabel('Current density (a.u.)')
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 1);
set(gca, 'linewidth', 1, 'FontSize', 16)

yyaxis right
% plot(Potential(1:801), 1000*CHI(800:end, 3)) % 20181116_D3_CHI
% plot(Potential(1:801), 1000*CHI(800:end, 2)) % Modified Au
plot(Potential(1:801), 1000*CHI(800:end, 1)) % Bare Au
ylim([-1.5 1.5])
ylabel('Current density (x 10^-^3 A/cm^2)')
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 1);
set(gca, 'linewidth', 1, 'FontSize', 16)
%% C1vsD3_drifted_subtract_1
x = (0:0.01:16)';  % 这里手动在potential的最后一行加上0，使其长度由1600变为1601。
xq = (0:0.005:16)'; 
figure('color', 'w');
potentialq1 = interp1(x,potential,xq);
plot(x, potential,'o', xq, potentialq1,':.'); % for examination
ylabel('Potential (V vs. Ag/AgCl)');
xlabel('Time (s)');

figure('color', 'w');

yyaxis left
plot(t, y3)
xlim([0 16]); ylim([-100 50])
xlabel('Time (s)')
ylabel('I_F - I_n_o_n_-_F (a.u.)')
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 1)
set(gca, 'linewidth', 1, 'FontSize', 16)

yyaxis right
plot(t, potentialq1)
ylim([-0.6 0]); ylabel('Potential (V vs. Ag/AgCl)')

%% for SI, the subtract process presented in pictures
load('E:\MoS2\20181116_MoS2_CH18-Au\TIFF_CV\Cut_D3-C1_mat.mat')
exp = cell(81, 1);
for ii = 1:81
    exp{ii, 1} = D3_C1.sampleArea{ii, 1}; % D3.sampleArea and C1.sampleArea
end

% figure('color', 'w');

savepath = 'E:\MoS2\20181116_MoS2_CH18-Au\TIFF_CV\Cut_D3_C1\';

for ii = 1:81
    % ii = 40;
    tif_ii = exp{ii, 1};
    
    if ii > 40
        Voltage = -0.4 + ((ii-41)/10)*0.1;
    else
        Voltage = -0.4 - ((ii-41)/10)*0.1;
    end
    
    localSums = imboxfilt(tif_ii, 13);
    localSums(localSums > 0) = 0;
    localSums(localSums < -1000) = -1000;
    
    imshow(localSums, 'DisplayRange',[], 'InitialMagnification', 'fit');
    title([num2str(Voltage), ' V'], 'FontSize', 14, 'FontWeight', 'bold');
    c = fire;
    c = flipud(c);
    map = colormap(c);
    % impixelinfo
    
    set(gca, 'CLim', [-1000 0]);
    h=colorbar;
    set(get(h,'title'),'string','Intensity (a.u.)', 'FontSize', 12);
    
    pause(0.05);
    saveas(gcf,[savepath, 'D3_C1_' num2str(ii, '%04d'), '.tif']);
end

%% Crop the .tif pictures
% Icropped = imcrop(I,rect) % rect = [197.5100   60.5100  655.9800  655.9800]

% imagesFolder = 'E:\MoS2\20181116_MoS2_CH18-Au\TIFF_CV\Cut_C1_fire';
% imagesFolder = 'E:\MoS2\20181116_MoS2_CH18-Au\TIFF_CV\Cut_D3_fire';
imagesFolder = 'E:\MoS2\20181116_MoS2_CH18-Au\TIFF_CV\Cut_D3_C1';

filePattern = [imagesFolder, '\*.tif'];
tifFiles = dir(filePattern);

% baseFileName = tifFiles(1).name;
% fullFileName = fullfile(imagesFolder, baseFileName);
% OriginalImage = imread(fullFileName);

% % Cropping Images and get the 'rect'
% [~, rect] = imcrop(OriginalImage);
rect = [197.5100   60.5100  655.9800  655.9800]; 
savepath = 'E:\MoS2\20181116_MoS2_CH18-Au\TIFF_CV\Cut_figureSI\';

for ii = 1:length(tifFiles)
    OriginalImage = imread(fullfile(imagesFolder, tifFiles(ii).name));
    croppedIMG = imcrop(OriginalImage, rect);
    filename = [savepath, 'D3_C1_' num2str(ii, '%04d'), '.tif'];
    imwrite(croppedIMG, filename);
end
%%

ii = 12;
tif_ii = exp{ii, 1};
% tif_ii = exp{ii, 1} - exp{9, 1}; % It is an optional.

Voltage = -0.3;

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
h = colorbar;
set(h, 'FontSize', 14, 'FontWeight', 'bold');
set(get(h,'title'),'string','[Ru(NH_3)_6]Cl_3^2^+ (mM)', 'FontSize', 12);
    


