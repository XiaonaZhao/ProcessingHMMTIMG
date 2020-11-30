%% 
exp = cell(11, 8);
fields = {'expName', 'tifPath', 'sROI', 'begin', 'pH', 'sampleRate', 'saveRoute', 'data'}; % the col number of exp
expTab = cell2struct(exp, fields, 2);

%%
prefix = ('H:\EDL\MoS2_20201021\_Timer\');
d = sortObj(dir([prefix, '*.mat']));
Fs = 100;
for ii = 1:size(expTab, 1)
    varMat = load([prefix, d(ii).name]);
    % varMat = load('G:\EDL\EDL_BareAu_20200813\_Timer\A1_data.mat');
    
    begin = triggerTime_PEIM(varMat.data, varMat.t, Fs);
    expTab(ii).begin = begin;
    expTab(ii).data = varMat.data;
    disp([prefix, d(ii).name]);
end

%%
prefix = ('E:\nona\Graphene_20201010_100x\_Result\sROI\');
d = sortObj(dir([prefix, '*.zip']));
for ii = 1:size(expTab,1)
    [expTab(ii).sROI] = ReadImageJROI([prefix, d(ii).name]);
    disp([prefix, d(ii).name]);
end

%%
load('F:\ZetaPotential\Graphene_20201010_100x\matlab_L1.mat')

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');
for m = 1:size(expTab, 1) 
% m = 8;
% for m = 1:6
    tic
    
    expName = expTab(m).expName;
    tifPath = expTab(m).tifPath;
    sROI = expTab(m).sROI;
    % mask = expTab(m).roiMask;
    begin = expTab(m).begin;
    saveRoute = expTab(m).saveRoute;
    % Fs = 100;
    data = expTab(m).data;
    voltage = data(:, 2);
    
%     [expTab(m).tif_Amp, expTab(m).tif_Ph, expTab(m).Amp, expTab(m).Phase] = MapGraphene(expName, tifPath, sROI, begin, saveRoute, voltage);
    [expTab(m).tif_Amp, expTab(m).tif_Ph, Amp, ~] = MapGraphene(expName, tifPath, sROI, begin, saveRoute, voltage);
    expTab(m).Amp_Au = Amp(1);
    expTab(m).Amp_Gra = Amp(2);
    
    processBar(size(expTab, 1), m, hwait)
    
    toc
end
delete(hwait);

cellpath = [saveRoute '\expTab_1.mat'];
save(cellpath, 'expTab', '-v7.3');
%%
% the charge density of the diffuse layer
ii = 1;
ref = 7;
Charge = (expTab(ii).tif_Amp).*(expTab(ii+4).tif_Amp)./(expTab(ii+4).tif_Amp - expTab(ii).tif_Amp);
% Charge = -imboxfilt(expTab(ii).tif_Amp, ref).*imboxfilt(expTab(ii+4).tif_Amp, ref)./(imboxfilt(expTab(ii).tif_Amp, ref) - imboxfilt(expTab(ii+4).tif_Amp, ref));
% Charge = imboxfilt(expTab(ii).tif_Amp, ref) - imboxfilt(expTab(ii+4).tif_Amp, ref);
% Charge = expTab(ii).tif_Amp - expTab(ii+4).tif_Amp;

% img = figure('color', 'w'); % amplitude
imshow(Charge, 'DisplayRange', [], 'InitialMagnification', 'fit');
title(expTab(ii).anion);
colormap default
h = colorbar;
set(get(h,'title'),'string','Amplitude (a.u.)');
% set(gca, 'CLim', [-10 10]);

%%
% figure('color', 'w'); % amplitude

ii = 8;
imshow(expTab(ii).tif_Amp, 'DisplayRange', [], 'InitialMagnification', 'fit');
title(expTab(ii).anion);
colormap default
h = colorbar;
% set(get(h,'title'),'string','Phase (radian)');
set(get(h,'title'),'string','Amplitude (a.u.)');
set(gca, 'CLim', [0 20]);
%%


