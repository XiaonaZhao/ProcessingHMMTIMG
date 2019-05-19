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
% TaS2_DC(expName, tifPath, maskPath, saveRoute);
% if CV
TaS2(expName, maskPath, loader, rate, saveRoute);

%%
load('F:\TaS2\20190324_TaS2_0318_ITO\result\expTab_TaS2_0324.mat')

for m = 1:size(expTab_TaS2, 1)
    expName = expTab_TaS2{m, 2};
    tifPath = expTab_TaS2{m, 3};
    maskPath = expTab_TaS2{m, 4};
    saveRoute = expTab_TaS2{m, 5};
    TaS2_DC(expName, tifPath, maskPath, saveRoute);
end

%% dc
[~, C2.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
C2.tifDir = dir(fullfile(C2.tifFile, '*.tiff'));
tif0 = double(imread(fullfile(C2.tifFile, C2.tifDir(1).name)));

[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
m =1;
[row, col] = ImageJroiLocation(sROI{m});

mask = ~imread('F:\TaS2\20190506_TaS2_ITO\Mask\C2.tif');
mask = mask((row(1):row(2)), (col(1):col(2)));

C2.beginFrame1 = 3082;
C2.endFrame1 = 3132;
C2.part1 = cell(51, 1);
for ii = C2.beginFrame1:C2.endFrame1
    tif  = double(imread(fullfile(C2.tifFile, C2.tifDir(ii).name))) - tif0;
    C2.part1{(ii-C2.beginFrame1) + 1, 1} = ii;
    C2.part1{(ii-C2.beginFrame1) + 1, 2} = tif((row(1):row(2)), (col(1):col(2)));
    C2.part1{(ii-C2.beginFrame1) + 1, 3} = mask.*C2.part1{(ii-C2.beginFrame1) + 1, 2};
end

for ii = 1:50
    C2.part2{ii, 1} = C2.part1{ii+1, 3} - C2.part1{ii, 3};
end

a = C2.part2{1, 1};
for ii = 1:10
    a = a + C2.part2{ii+1, 1};
end
a = a/11;
b = length(find(a(:)~=0));
c =  sum(a(:))/b;

%%
ii = 26;
d = C2.part1{ii, 3};
d = imboxfilt(d, 3);
d(d > 0) = 0;
d = d.*mask;
imshow(d, 'DisplayRange',[], 'InitialMagnification', 'fit');
% colormap jet 
% impixelinfo

%%

C2.part1diff = cell(size(C2.part1, 1) -1, 1);
for ii = 1:length(C2.part1diff)
    C2.part1diff{ii} = C2.part1{ii+1, 2} - C2.part1{ii, 2};
end


temp = C2.part1diff;
C2.part1_3D = zeros([size(img) 349]);
for n = 1:length(temp)
    C2.part1_3D(:, :, n) = temp{n, 1};
end
clear temp

[C2.part1diff_value, C2.part1diff_time] = imgSeqDiff_s(C2.part1_3D);

%% 20190501_TaS2_ITO
load('H:\TaS2\20190507_TiS2_ITO\Result\Darker\Darker.mat')

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');
% for m = 1:size(expTab, 1)
for m = [1 2 3 4 6 7]
    expName = expTab(m).expName;
    tifPath = expTab(m).tifPath;
    Mask = expTab(m).Mask;
    begin = expTab(m).begin;
    saveRoute = expTab(m).saveRoute;
    rate = expTab(m).ScanRate;
    
    TaS2_batch(expName, tifPath, Mask, begin, rate, saveRoute);
    
    disp(['The latest progress is about ' num2str(m) '.']);
    processBar(size(expTab, 1), m, hwait)
    
end

%% 20190501_TaS2_ITO
% 100 mV/s
load('H:\TaS2\20190507_TiS2_ITO\Result\Lighter\B1.mat')
expTab(1).potential = Value.potential(1601:3200, 1);
expTab(1).dIntensity = Value.outside(1601:3200, 1:2);
expTab(1).dROImean = Value.dROImean(1601:3200, 1);
% 200 mV/s
load('H:\TaS2\20190507_TiS2_ITO\Result\Lighter\B2.mat')
expTab(2).potential = Value.potential(801:1600, 1);
expTab(2).dIntensity = Value.outside(801:1600, 1:2);
expTab(2).dROImean = Value.dROImean(801:1600, 1);
% 300 mV/s
load('H:\TaS2\20190507_TiS2_ITO\Result\Lighter\B3.mat')
expTab(3).potential = Value.potential(536:1069, 1);
expTab(3).dIntensity = Value.outside(536:1069, 1:2);
expTab(3).dROImean = Value.dROImean(536:1069, 1);
% 400 mV/s
load('H:\TaS2\20190507_TiS2_ITO\Result\Lighter\B4.mat')
expTab(4).potential = Value.potential(401:800, 1);
expTab(4).dIntensity = Value.outside(401:800, 1:2);
expTab(4).dROImean = Value.dROImean(401:800, 1);
% 500 mV/s
load('H:\TaS2\20190507_TiS2_ITO\Result\Lighter\B5.mat')
expTab(5).potential = Value.potential(321:640, 1);
expTab(5).dIntensity = Value.outside(321:640, 1:2);
expTab(5).dROImean = Value.dROImean(321:640, 1);
% 50 mV/s
load('H:\TaS2\20190507_TiS2_ITO\Result\Lighter\B6.mat')
expTab(6).potential = Value.potential(3201:6400, 1);
expTab(6).dIntensity = Value.outside(3201:6400, 1:2);
expTab(6).dROImean = Value.dROImean(3201:6400, 1);
% % 50 mV/s
% load('I:\TaS2\20190501_TaS2_ITO\Result\B6_2.mat')
% expTab(7).potential = Value.potential(3201:6400, 1);
% expTab(7).dIntensity = Value.outside(3201:6400, 1:2);
% expTab(7).dROImean = Value.dROImean(3201:6400, 1);
% % 100 mV/s
% load('I:\TaS2\20190501_TaS2_ITO\Result\B7.mat')
% expTab(8).potential = Value.potential(1601:3200, 1);
% expTab(8).dIntensity = Value.outside(1601:3200, 1:2);
% expTab(8).dROImean = Value.dROImean(1601:3200, 1);
% % 100 mV/s
% load('I:\TaS2\20190501_TaS2_ITO\Result\B8.mat')
% expTab(9).potential = Value.potential(1601:3200, 1);
% expTab(9).dIntensity = Value.outside(1601:3200, 1:2);
% expTab(9).dROImean = Value.dROImean(1601:3200, 1);
%%
figure('color', 'w');
hold on
for ii = 1:size(expTab, 1)
    plot(expTab(ii).num, expTab(ii).ROImean)
end
hold off
%%
figure('color', 'w');
hold on
plot(expTab(6).potential, 0.8*expTab(6).dIntensity, 'color', [0.6350, 0.0780, 0.1840]) % 50 mV/s
plot(expTab(1).potential, expTab(1).dIntensity, 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).potential, expTab(2).dIntensity, 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
plot(expTab(3).potential, expTab(3).dIntensity, 'color', [0.4660, 0.6740, 0.1880]) % 300 mV/s
plot(expTab(4).potential, expTab(4).dIntensity, 'color', [0.3010, 0.7450, 0.9330]) % 400 mV/s
plot(expTab(5).potential, expTab(5).dIntensity, 'color', [0, 0.4470, 0.7410]) % 500 mV/s
xlabel('Potential (V)');
ylabel('\DeltaIntensity'' (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%%
figure('color', 'w');
hold on
plot(expTab(6).potential, expTab(6).dROImean, 'color', [0.6350, 0.0780, 0.1840]) % 50 mV/s
plot(expTab(1).potential, 2*expTab(8).dROImean, 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).potential, expTab(2).dROImean, 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
plot(expTab(3).potential, expTab(3).dROImean, 'color', [0.4660, 0.6740, 0.1880]) % 300 mV/s
plot(expTab(4).potential, expTab(4).dROImean, 'color', [0.3010, 0.7450, 0.9330]) % 400 mV/s
plot(expTab(5).potential, expTab(5).dROImean, 'color', [0, 0.4470, 0.7410]) % 500 mV/s
xlabel('Potential (V)');
ylabel('Averaged \DeltaIntensity'' (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%%
figure('color', 'w');
hold on
plot(expTab(6).potential, 0.52*expTab(6).ROImean(1:3200,1), 'color', [0.6350, 0.0780, 0.1840]) % 50 mV/s
plot(expTab(1).potential, expTab(1).ROImean(1:1600,1), 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).potential, expTab(2).ROImean(1:800,1), 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
plot(expTab(3).potential, expTab(3).ROImean(1:534,1), 'color', [0.4660, 0.6740, 0.1880]) % 300 mV/s
plot(expTab(4).potential, expTab(4).ROImean(1:400,1), 'color', [0.3010, 0.7450, 0.9330]) % 400 mV/s
plot(expTab(5).potential, expTab(5).ROImean(1:320,1), 'color', [0, 0.4470, 0.7410]) % 500 mV/s
xlabel('Potential (V)');
ylabel('Intensity (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%%
figure('color', 'w');
hold on
plot(expTab(6).num, -0.52*expTab(6).ROImean, 'color', [0.6350, 0.0780, 0.1840]) % 50 mV/s
plot(expTab(1).num, -expTab(1).ROImean, 'color', [0.8500, 0.3250, 0.0980]) % 100 mV/s
plot(expTab(2).num, -expTab(2).ROImean, 'color', [0.9290, 0.6940, 0.1250]) % 200 mV/s
plot(expTab(3).num, -expTab(3).ROImean, 'color', [0.4660, 0.6740, 0.1880]) % 300 mV/s
plot(expTab(4).num, -expTab(4).ROImean, 'color', [0.3010, 0.7450, 0.9330]) % 400 mV/s
plot(expTab(5).num, -expTab(5).ROImean, 'color', [0, 0.4470, 0.7410]) % 500 mV/s
xlabel('Potential (V)');
ylabel('Intensity (a.u.)');
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off
%% Test the
% figure('color', 'w');
curve = lowp(Value.ROImean, 5, 22, 0.1, 20, 100); % Bright Field, CV;
plot(X, diff(Value.ROImean),  X, diff(curve))
% set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
% xlim([2300 2700])
% ylim([100 450])

%%
varMat = load('H:\TaS2\20190507_TiS2_ITO\Timer\A1_data.mat');
begin = triggerTime(varMat.data, varMat.t);
expTab(1).begin = begin;

%% 20190506_TaS2_ITO

prefix = ('F:\TaS2\20190506_TaS2_ITO\Result\');
d = dir([prefix, '*.mat']);
for ii = 1:size(expTab, 1)
    load([prefix, d(ii).name]) ;
    expTab(ii).potential = Value.potential(1601:3200, 1);
    expTab(ii).dIntensity = Value.outside(1601:3200, 1:2);
    expTab(ii).dROImean = Value.dROImean(1:1600, 1);
end

figure('color', 'w');
hold on
plot(expTab(2).potential, expTab(2).dIntensity, 'color', [0, 0.4470, 0.7410]) % K
plot(expTab(4).potential, expTab(4).dIntensity, 'color', [0.3010, 0.7450, 0.9330]) % Na
plot(expTab(6).potential, expTab(6).dIntensity, 'color', [0.4660, 0.6740, 0.1880]) % Li
plot(expTab(8).potential, expTab(8).dIntensity, 'color', [0.9290, 0.6940, 0.1250]) % Mg
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off

figure('color', 'w');
hold on
plot(expTab(2).potential, expTab(2).dROImean, 'color', [0, 0.4470, 0.7410]) % K 3.31
plot(expTab(4).potential, expTab(4).dROImean, 'color', [0.3010, 0.7450, 0.9330]) % Na 3.58
plot(expTab(6).potential, expTab(6).dROImean, 'color', [0.4660, 0.6740, 0.1880]) % Li 3.82
plot(expTab(8).potential, expTab(8).dROImean, 'color', [0.9290, 0.6940, 0.1250]) % Mg 4.28
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
hold off

%% Rename image Seq
folder_name = uigetdir;
folder = dir(folder_name);

oldname = cell(length(folder)-2, 1);
for ii = 3:length(folder)
   oldname{ii-2} = folder(ii).name;
end

newname = cell(length(oldname), 1);
for ii = 2:length(oldname)
   a = oldname{ii};
   b = str2num(a(12:end-6));
   c = num2str(b, '%06d');
   newname{ii} = ['B4_', c, '.tiff'];
   movefile([folder_name '\' oldname{ii}], [folder_name '\' newname{ii}])
end

%%
saveRoute = 'F:\TaS2\20190506_TaS2_ITO\Result\Pic';
% V = [-0.72 -0.701 -0.682 -0.64 -0.604 -0.575]';
V = [-0.686 -0.675 -0.631]';
n = zeros(size(expTab, 1), length(V));
m = zeros(size(expTab, 1), length(V));
for ii = [2 3 6]
    tifDir = dir(fullfile(expTab(ii).tifPath, '*.tiff'));
    tif0 = double(imread(fullfile(expTab(ii).tifPath, tifDir(1).name)));

    a = expTab(ii).potential;
    for jj = 1:length(V)
        b = a(:) - V(jj);
        [~, n(ii, jj)] = min(abs(b));
        m(ii, jj) = n(ii, jj) + expTab(ii).begin.frame +length(expTab(ii).potential);
        tif = double(imread(fullfile(expTab(ii).tifPath, tifDir(m(ii, jj)).name))) - tif0;
        tif1 = tif((row(1):row(2)), (col(1):col(2)));
        figPath = [saveRoute '\B' num2str(ii) '_' num2str(jj) '.tif'];
        imwrite(tif1, figPath)
    end
end
%%

