function [Amp, Pha] = MoS2_charging(expName, tifPath, sROI, begin, saveRoute, voltage)

L = 1024;
Fs = 100;
voltage = voltage(round(begin.CS1/100):(round(begin.CS1/100) + L*Fs - 1));

Value.tifFile = tifPath;
Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));
Value.validDir = Value.tifDir(begin.frame:(begin.frame + L -1));

[row, col] = ImageJroiLocation(sROI{1});
%%
tif0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
tif = cell(L, 1); % why not 'tif = cell(L, 1);'
for ii = 1:L
    temp = (double(imread(fullfile(Value.tifFile, Value.validDir(ii).name))) - tif0)./tif0;
    tif{ii, 1} = temp(row(1):row(2), col(1):col(2));
end

tif_Amp_peak = zeros(size(tif{1}));
tif_Ph = zeros(size(tif{1}));

tif_3D = zeros(size(tif{1}, 1), size(tif{1}, 2), L);
parfor ii = 1:L
    tif_3D(:, :, ii) = tif{ii, 1};
end

for ii = 1:size(tif{1}, 1)
    parfor jj = 1:size(tif{1}, 2)
        X = reshape(tif_3D(ii, jj, :), L, 1);
        Y = fft(X);
        P2 = abs(Y/L);
        tif_Amp_peak(ii, jj) = max(2*P2(2:end-1));
        [~, tif_Ph(ii, jj)] = freqrespmeasure(voltage, X);
    end
end
% Value.A_peak = tif_Amp_peak;
% Value.tif_Ph = tif_Ph;
% clear tif_3D

img = figure('color', 'w');
imshow(tif_Amp_peak, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap default
h=colorbar;
set(get(h,'title'),'string','Intensity (a.u.)', 'FontSize', 14);
set(gca, 'CLim', [0 0.05]);
%%
figPath = [saveRoute '\' expName];
saveas(img, figPath, 'fig')

% mask = ~imread(mask);
% Value.Amplitude = ROImean(tif_Amp_peak, mask);
% Amp = Value.Amplitude;

% intensity = zeros(L, 1);
% for ii = 1:L
%     intensity(ii ,1) = ROImean(tif{ii, 1}, mask);
% end
% Value.intensity = intensity;
% clear tif
% [Amp, Pha] = freqrespmeasure(voltage, intensity);
% 
% cellpath = [saveRoute '\' expName '.mat']; 
% save(cellpath, 'Value', '-v7.3');

end