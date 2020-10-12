% for Graphene_20200108
tic

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');

load('G:\EDL\Graphene_20200108\_Result\matlab.mat')
Amp = zeros(size(expTab, 1), 2);
for m = 1:size(expTab, 1)
    expName = expTab(m).expName;
    tifPath = expTab(m).tifPath;
    begin = expTab(m).begin;
    sROI = expTab(m).sROI;
    [Amp(m, 1), Amp(m, 2)] = edl_1(tifPath, sROI, begin);
    
    DeltaV = 0.04; % 10 mVpp
    Ca = Amp/DeltaV;
    
    disp([expName ': C1 = ' num2str(Ca(m, 1)) ', C2 = ' num2str(Ca(m, 2))]);
    processBar(size(expTab, 1), m, hwait)
end

delete(hwait);

toc



function [Amp1, Amp2] = edl_1(tifPath, sROI, begin)
L = 1024;

tifDir = dir(fullfile(tifPath, '*.tiff'));
validDir = tifDir(begin.frame:(begin.frame + 2*L -1));
% Value.validDir = Value.tifDir(begin.frame + L:(begin.frame + 2*L -1)); %
% for the second part

[row, col] = ImageJroiLocation(sROI);

tifValue = zeros(2*L, 1);
tif0 = double(imread(fullfile(tifPath, tifDir(1).name)));
for ii = 1:2*L
    % temp = (double(imread(fullfile(tifPath, validDir(ii).name))) - tif0)./tif0;
    temp = double(imread(fullfile(tifPath, validDir(ii).name))) - tif0;
    temp = temp((row(1):row(2)), (col(1):col(2)));
    tifValue(ii, 1) = mean(temp(:));
end
Y = fft(tifValue(1:L));
P2 = abs(Y/L);
Amp1 = max(2*P2(4:ceil(L/2+1)-1));

Y = fft(tifValue(L+1:2*L));
P2 = abs(Y/L);
Amp2 = max(2*P2(4:ceil(L/2+1)-1));

end