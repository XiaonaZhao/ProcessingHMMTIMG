% for Graphene_20200827
%About the ROIs
%   Mask_intensity is for the Mask image.
%   MoS2 is the rectangle ROI of the sample.
%   MoS2 is the rectangle ROI of the sample.

tic

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');

load('G:\EDL\MoS2_20200827\_Result\matlab_B.mat')
L_i = zeros(11, 4, 3);
for mm = 1:size(expTab, 1)
    expName = expTab(mm).expName;
    tifPath = expTab(mm).tifPath;
    Mask = expTab(mm).Mask;
    begin = expTab(mm).begin;
    sROI = expTab(mm).sROI;
    [Mask_i, MoS2, Au, small] = edl_1(tifPath, Mask, sROI, begin);
    
    DeltaV = 0.1; % 25 mVpp
    for nn = 1:3
        L_i(mm, 1, nn) = Mask_i(nn)/DeltaV;
        L_i(mm, 2, nn) = MoS2(nn)/DeltaV;
        L_i(mm, 3, nn) = Au(nn)/DeltaV;
        L_i(mm, 4, nn) = small(nn)/DeltaV;
    end
    expTab(mm).Mask_intensity = Mask_i;
    expTab(mm).MoS2 = MoS2;
    expTab(mm).Au = Au;
    expTab(mm).small = small;
    
    disp([expName ' finished.']);
    processBar(size(expTab, 1), mm, hwait)
end

delete(hwait);

toc


function [Mask_i, MoS2, Au, small] = edl_1(tifPath, Mask, sROI, begin)
L = 1024;

tifDir = dir(fullfile(tifPath, '*.tiff'));
validDir = tifDir(begin.frame:(begin.frame + 3*L -1));
tif0 = double(imread(fullfile(tifPath, tifDir(1).name)));

Mask = imread(Mask);
Mask = ~Mask;
intensity_M = zeros(3*L, 1);
for ii = 1:3*L
    % temp = (double(imread(fullfile(tifPath, validDir(ii).name))) - tif0)./tif0.*Mask;
    temp = (double(imread(fullfile(tifPath, validDir(ii).name))) - tif0).*Mask;
    intensity_M(ii, 1) = sum(temp(:))/sum(Mask(:));
end
P2 = abs(fft(intensity_M(1:L))/L);
Mask_i(1) = max(2*P2(4:ceil(L/2+1)-1));
P2 = abs(fft(intensity_M(L+1:2*L))/L);
Mask_i(2) = max(2*P2(4:ceil(L/2+1)-1));
P2 = abs(fft(intensity_M(2*L+1:3*L))/L);
Mask_i(3) = max(2*P2(4:ceil(L/2+1)-1));

sROI_i = zeros(3, 3);
for ii = 1:3 % sROI for MoS2, Au, small
    [row, col] = ImageJroiLocation(sROI{ii});
    
    tifValue = zeros(3*L, 1);
    for jj = 1:3*L % L, 2L, 3L
        % temp = (double(imread(fullfile(tifPath, validDir(jj).name))) - tif0)./tif0;
        temp = double(imread(fullfile(tifPath, validDir(jj).name))) - tif0;
        temp = temp((row(1):row(2)), (col(1):col(2)));
        tifValue(jj, 1) = mean(temp(:));
    end
    P2 = abs(fft(tifValue(1:L))/L);
    sROI_i(1, ii) = max(2*P2(4:ceil(L/2+1)-1));
    P2 = abs(fft(tifValue(L+1:2*L))/L);
    sROI_i(2, ii) = max(2*P2(4:ceil(L/2+1)-1));
    P2 = abs(fft(tifValue(2*L+1:3*L))/L);
    sROI_i(3, ii) = max(2*P2(4:ceil(L/2+1)-1));
    
end
MoS2 = sROI_i(1, :);
Au = sROI_i(2, :);
small = sROI_i(3, :);

end