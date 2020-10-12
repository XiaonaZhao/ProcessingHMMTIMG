% forMoS20191121

tic

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');

% load('G:\EDL\MoS2_20191121\_Result\matlab_C.mat') % Load C
% load('G:\EDL\MoS2_20191121\_Result\matlab_D.mat') % Load D
load('G:\EDL\MoS2_20191121\_Result\matlab_E.mat') % Load E

for mm = 1:size(expTab, 1)
    expName = expTab(mm).expName;
    tifPath = expTab(mm).tifPath;
    begin = expTab(mm).begin;
    
    [MoS2, Au] = edl_1(tifPath, sROI, begin);
    
    DeltaV = 0.4; % 100 mVpp
    expTab(mm).MoS2 = MoS2/DeltaV;
    expTab(mm).Au = Au/DeltaV;
    
    disp([expName ' finished.']);
    processBar(size(expTab, 1), mm, hwait)
end

delete(hwait);

toc


function [MoS2, Au] = edl_1(tifPath, sROI, begin)
L = 1024;

tifDir = dir(fullfile(tifPath, '*.tiff'));
validDir = tifDir(begin.frame:(begin.frame + L -1));
tif0 = double(imread(fullfile(tifPath, tifDir(1).name)));

tifValue = zeros(L, length(sROI));
for ii = 1:L % L
    % tif = (double(imread(fullfile(tifPath, validDir(ii).name))) - tif0)./tif0;
    tif = double(imread(fullfile(tifPath, validDir(ii).name))) - tif0;
    for jj = 1:length(sROI) % sROI for Gra1, Gra2, Au
        [row, col] = ImageJroiLocation(sROI{jj});
        temp = tif((row(1):row(2)), (col(1):col(2)));
        tifValue(ii, jj) = mean(temp(:));
    end
end

sROI_i = zeros(1, length(sROI));
for jj = 1:length(sROI)
    P2 = abs(fft(tifValue(1:L, jj))/L);
    sROI_i(1, jj) = max(2*P2(3:ceil(L/2+1)-1));
end

MoS2 = sROI_i(1);
Au = sROI_i(2);

end