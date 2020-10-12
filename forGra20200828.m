% forGra20200828

tic

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');

load('G:\EDL\Graphene_20200828\_Result\matlab_B.mat') % Load expTab, sROI

L_i = zeros(11, 3, 3);
for mm = 1:size(expTab, 1)
    expName = expTab(mm).expName;
    tifPath = expTab(mm).tifPath;
    begin = expTab(mm).begin;
    
    [Gra1, Gra2, Au] = edl_1(tifPath, sROI, begin);
    
    DeltaV = 0.1; % 25 mVpp
    for nn = 1:3
        L_i(mm, 1, nn) = Gra1(nn)/DeltaV;
        L_i(mm, 2, nn) = Gra2(nn)/DeltaV;
        L_i(mm, 3, nn) = Au(nn)/DeltaV;
    end
    expTab(mm).Gra1 = Gra1;
    expTab(mm).Gra2 = Gra2;
    expTab(mm).Au = Au;
    
    disp([expName ' finished.']);
    processBar(size(expTab, 1), mm, hwait)
end

delete(hwait);

toc


function [Gra1, Gra2, Au] = edl_1(tifPath, sROI, begin)
L = 1024;

tifDir = dir(fullfile(tifPath, '*.tiff'));
validDir = tifDir(begin.frame:(begin.frame + 3*L -1));
tif0 = double(imread(fullfile(tifPath, tifDir(1).name)));

tifValue = zeros(3*L, length(sROI));
for ii = 1:3*L % L, 2L, 3L
    % tif = (double(imread(fullfile(tifPath, validDir(ii).name))) - tif0)./tif0;
    tif = double(imread(fullfile(tifPath, validDir(ii).name))) - tif0;
    for jj = 1:length(sROI) % sROI for Gra1, Gra2, Au
        [row, col] = ImageJroiLocation(sROI{jj});
        temp = tif((row(1):row(2)), (col(1):col(2)));
        tifValue(ii, jj) = mean(temp(:));
    end
end

sROI_i = zeros(3, length(sROI));
for jj = 1:length(sROI)
    P2 = abs(fft(tifValue(1:L, jj))/L);
    sROI_i(1, jj) = max(2*P2(4:ceil(L/2+1)-1));
    P2 = abs(fft(tifValue(L+1:2*L, jj))/L);
    sROI_i(2, jj) = max(2*P2(4:ceil(L/2+1)-1));
    P2 = abs(fft(tifValue(2*L+1:3*L, jj))/L);
    sROI_i(3, jj) = max(2*P2(4:ceil(L/2+1)-1));
end

Gra1 = sROI_i(1, :);
Gra2 = sROI_i(2, :);
Au = sROI_i(3, :);

end