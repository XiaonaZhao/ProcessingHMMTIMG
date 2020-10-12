function [tif_Amp, tif_Ph, Amp, Phase] = MapGraphene(expName, tifPath, sROI, begin, saveRoute, voltage)

L = 1024;
Fs = 100; % sample rate
voltage = voltage(begin.CS1:(begin.CS1 + L*Fs - 1)); % changed
% voltage = voltage((begin.CS1 + L*10000/Fs):(begin.CS1 + 3*L*10000/Fs - 1));

tifDir = dir(fullfile(tifPath, '*.tiff'));
validDir = tifDir(begin.frame:(begin.frame + L -1));  % changed
% validDir = tifDir((begin.frame + L):(begin.frame + 2*L -1));
tif0 = double(imread(fullfile(tifPath, tifDir(1).name)));
tif = cell(L, 1); % why not 'tif = cell(L, 1);'
tifValue = zeros(L, length(sROI));

for ii = 1:L
    % tif{ii, 1} = (double(imread(fullfile(tifPath, validDir(ii).name))) - tif0)./tif0;
    tif{ii, 1} = double(imread(fullfile(tifPath, validDir(ii).name))) - tif0;
    temp = tif{ii, 1};
    for jj = 1:length(sROI)
        [row, col] = ImageJroiLocation(sROI{jj});
        temp2 = temp(row(1):row(2), col(1):col(2));
        tifValue(ii, jj) = mean(temp2(:));
    end
    
end

tif_3D = zeros(size(tif{1}, 1), size(tif{1}, 2), L);
parfor ii = 1:L
    tif_3D(:, :, ii) = tif{ii, 1};
end

tif_Amp = zeros(size(tif{1}));
tif_Ph = zeros(size(tif{1}));

for ii = 1:size(tif{1}, 1)
    parfor jj = 1:size(tif{1}, 2)
        X = reshape(tif_3D(ii, jj, :), L, 1);
        Y = fft(X);
        P2 = abs(Y/L);
        tif_Amp(ii, jj) = max(2*P2(3:end-1));
        [~, tif_Ph(ii, jj)] = freqrespmeasure(voltage, X);
    end
    
end

img = figure('color', 'w'); % amplitude
imshow(tif_Amp, 'DisplayRange', [], 'InitialMagnification', 'fit');
title(expName);
colormap default
h = colorbar;
set(get(h,'title'),'string','Amplitude (a.u.)');
% set(gca, 'CLim', [0 0.05]);
figPath = [saveRoute '\' expName '_amp'];
saveas(img, figPath, 'fig')
close all

img = figure('color', 'w'); % phase
imshow(tif_Ph, 'DisplayRange', [], 'InitialMagnification', 'fit');
title(expName);
colormap default
h = colorbar;
set(get(h,'title'),'string','Phase (radian)');
figPath = [saveRoute '\' expName '_ph'];
saveas(img, figPath, 'fig')
close all

Amp = zeros(length(sROI), 1);
Phase = zeros(length(sROI), 1);
for ii = 1:length(sROI)
    Y = fft(tifValue(:, ii));
    P2 = abs(Y/L);
    Amp(ii) = max(2*P2(3:end-1));
%     [~, Phase(ii)] = freqrespmeasure(voltage, tifValue(:, ii));
end

end