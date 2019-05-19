function TaS2_batch(expName, tifPath, Mask, begin, rate, saveRoute)
% ITO

tic

% All raw images
Value.tifFile = tifPath;
Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));

% Potential and ScanRate
if rate == 300
    r = -0.003;
    potential1 = (0 : r : -0.8)';
    potential2 = (-0.799 : (-r) : 0)';
    potential3 = (0.002 : r : -0.8)';
    potential4 = (-0.798 : (-r) : 0)';
    Value.potential = [potential1' potential2' potential3' potential4']';
    clear potential1 potential2 potential3 potential4
else
    
    if rate == 50
        r = -0.0005;
    elseif rate == 100
        r = -0.001;
    elseif rate == 200
        r = -0.002;
    elseif rate == 400
        r = -0.004;
    elseif rate == 500
        r = -0.005;
    else
        return
    end
    
    potential1 = (0 : r : -0.8)';
    potential2 = ((-0.8-r) : (-r) : -0.0000001)';
    Value.potential = [potential1' potential2' potential1' potential2']';
    clear potential1 potential2
end

% Valid image sequence
Value.validDir = Value.tifDir(begin.frame:(begin.frame+length(Value.potential)));
Value.begin = begin;

% for n = 1:length(Value.maskNames)
Mask = imread(Mask);
mask = ~Mask;

if sum(mask(:)) == 0
    return
end

points = ReadTifMaskPoint(Value.tifFile, Value.validDir, mask);

Fs = 100; % sampling rate
col = size(points, 2);
curve = zeros(size(points));
for ii = 1:1:col
    % curve(:, ii) = lowp(points(:, ii), 1, 36, 0.1, 20, Fs); % SPR, 20;
    curve(:, ii) = lowp(points(:, ii), 2, 11, 0.1, 20, Fs); % Bright Field, CV;
end
clear points

X = (1:1:(size(curve, 1)-1))';
dcurve = zeros(size(curve, 1)-1, size(curve, 2));

for ii = 1:col
    dcurve(:, ii) = diff(curve(:, ii));
end

Value.outside = -figSketch(dcurve);
outside = Value.outside;
img = figure('color','w');
hold on
for ii = 1:2
    plot(X, outside(:, ii), '.k')
end
xlabel('Frames'); ylabel('\DeltaIntensity''');
title([expName ' K_2SO_4 ROI, ' num2str(rate) ' mV/s'])
hold off
figPath = [saveRoute '\' expName '_roi'];
saveas(img, figPath, 'fig')
% end

img2 = figure('color','w');
hold on
% for n = 1:length(Value.maskNames)
% outside = Value.outside;
for ii = 1:2
    plot(Value.potential, outside(:, ii), '.k')
end
% end
xlabel('Potential/V'); ylabel('\DeltaIntensity''');
title([expName ' \DeltaIntensity'' with Potential, K_2SO_4, ' num2str(rate) ' mV/s'])
hold off
figPath2 = [saveRoute '\' expName '_intensityVSpotential'];
saveas(img2, figPath2, 'fig')


% % ROIMEAN
tif0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
for ii = Value.begin.frame:(Value.begin.frame+length(Value.potential))
    tif  = double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name))) - tif0;
%     Value.ROImean((ii-Value.begin.frame+1), 1) = ROImean(tif, mask); % TaS2
    Value.ROImean((ii-Value.begin.frame+1), 1) = -ROImean(tif, mask); % TiS2
end
temp = lowp(Value.ROImean, 2, 12, 0.1, 20, 100); % TaS2
% temp = lowp(Value.ROImean, 5, 22, 0.1, 20, 100); % TiS2 is not good
% enough
Value.dROImean = -diff(temp); clear temp
img3 = figure('color','w');
plot(Value.potential, Value.dROImean, 'k')
xlabel('Potential/V'); ylabel('\DeltaIntensity''');
title([expName ' Averagered Intensity'' with Potential, K_2SO_4, ' num2str(rate) ' mV/s'])
hold off
figPath3 = [saveRoute '\' expName '_AveragedintensityVSpotential'];
saveas(img3, figPath3, 'fig')


close all

% Save Value
cellpath = [saveRoute '\' expName '.mat'];
save(cellpath, 'Value', '-v7.3');

toc

end