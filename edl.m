function Ampl = edl(expName, tifPath, sROI, begin, saveRoute)
% This is a formula for calculate electrochemical double layer capacitance


L = 1024;

Value.tifFile = tifPath;
Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));
Value.validDir = Value.tifDir(begin.frame:(begin.frame + L -1));
% Value.validDir = Value.tifDir(begin.frame + L:(begin.frame + 2*L -1)); %
% for the second part


tif0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
for ii = begin.frame:(begin.frame + L -1)
    temp = double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name)));
    Value.tif{(ii-begin.frame+1), 1} = (temp - tif0)./tif0;
end
clear temp

tif_A_peak = zeros(size(tif0, 1), size(tif0, 2));

tif_3D = zeros(size(tif0, 1), size(tif0, 2), L);
for ii = 1:L
    tif_3D(:, :, ii) = Value.tif{ii, 1};
end

for ii = 1:size(tif0, 1)
    parfor jj = 1:size(tif0, 2)
        X = reshape(tif_3D(ii, jj, :), L, 1);
        Y = fft(X);
        P2 = abs(Y/L);
        tif_A_peak(ii, jj) = max(2*P2(2:end-1));
    end
end
Value.A_peak = tif_A_peak;
clear tif_3D

img = figure('color', 'w');
imshow(tif_A_peak, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet
figPath = [saveRoute '\' expName];
saveas(img, figPath, 'fig')
figName = [saveRoute '\' expName '.tif'];
imwrite(tif_A_peak, figName);

Ampl = zeros(length(sROI), 1);
maxAmpl = cell(length(sROI), 1);
for mm = 1:length(sROI)
    [row, col] = ImageJroiLocation(sROI{mm});
    
    tifValue = zeros(L, 1);
    for ii = begin.frame:(begin.frame + L -1)
        temp = Value.tif{(ii-begin.frame+1), 1};
        temp = temp((row(1):row(2)), (col(1):col(2)));
        tifValue(ii, 1) = mean(temp(:));
    end
    clear temp temp1
    Y = fft(tifValue);
    P2 = abs(Y/L);
    Ampl(mm, 1) = max(2*P2(2:end-1));
    
    temp = tif_A_peak((row(1):row(2)), (col(1):col(2)));
    maxAmpl{mm, 1} = temp(:);
    
end
Value.ampl = Ampl;
Value.maxAmpl = maxAmpl;

cellpath = [saveRoute '\' expName '.mat']; 
save(cellpath, 'Value', '-v7.3');

close all


end
