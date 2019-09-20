function A = MoS2_charging(expName, tifPath, mask, begin, saveRoute)

L = 1024;

Value.tifFile = tifPath;
Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));
Value.validDir = Value.tifDir(begin.frame:(begin.frame+length(Value.potential)));

tif0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
for ii = Value.begin.frame:(Value.begin.frame+L)
    Value.tif{ii, 1}  = (double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name))) - tif0)./tif0;
end

tif_A_peak = zeros(size(tif0, 1), size(tif0, 2));

tif_3D = zeros(480, 640, 1024);
for ii = 1:L
    tif_3D(:, :, ii) = Value.tif{ii, 1};
end

for ii = 1:size(tif0, 1)
    for jj = 1:size(tif0, 2)
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

mask = ~imread(mask);
Value.Amplitude = ROImean(tif_A_peak, mask);
A = Value.Amplitude;

cellpath = [saveRoute '\' expName '.mat']; 
save(cellpath, 'Value', '-v7.3');

end