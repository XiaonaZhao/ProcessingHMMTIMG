function TaS2(expName, maskPath, loader, rate, saveRoute)

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
    else
        return
    end
    
    potential1 = (0 : r : -0.8)';
    potential2 = ((-0.8-r) : (-r) : -0.0000001)';
    Value.potential = [potential1' potential2' potential1' potential2']';
    clear potential1 potential2
end

[~, Value.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));

varMat = load(loader);
begin = triggerTime(varMat.data, varMat.t);
Value.validDir = Value.tifDir(begin.frame:(begin.frame+length(Value.potential)));
Value.begin = begin;

Value.maskPath = maskPath;
[~, Value.maskNames] = ReadTifFileNames(Value.maskPath);

for n = 1:length(Value.maskNames)
    Mask = imread(fullfile(Value.maskPath, Value.maskNames{n}));
    mask = ~Mask;
    % span = 8;
    if sum(mask(:)) == 0
        return
    end
    
    points = ReadTifMaskPoint(Value.tifFile, Value.validDir, mask);
    
    Fs = 106;
    col = size(points, 2);
    curve = zeros(size(points));
    for ii = 1:1:col
        curve(:, ii) = lowp(points(:, ii), 1, 36, 0.1, 20, Fs); % SPR, 20;
        %         curve(:, ii) = lowp(points(:, ii), 2, 11, 0.1, 20, Fs); % Bright Field, CV;
    end
    clear points
    
    X = (1:1:(size(curve, 1)-1))';
    dcurve = zeros(size(curve, 1)-1, size(curve, 2));
    
    for ii = 1:col
        dcurve(:, ii) = diff(curve(:, ii));
    end
    
    Value.outside{n, 1} = figSketch(dcurve);
    outside = Value.outside{n, 1};
    img = figure('color','w');
    hold on
    for ii = 1:2
        plot(X, outside(:, ii), '.k')
    end
    xlabel('Frames'); ylabel('\DeltaIntensity''');
    title([expName ' Na_2SO_4 ROI' num2str(n)])
    hold off
    figPath = [saveRoute '\' expName '_roi' num2str(n) ];
    saveas(img, figPath, 'fig')
    
end

cellpath = [saveRoute '\' expName '.mat'];
save(cellpath, 'Value', '-v7.3');

img2 = figure('color','w');
hold on
for n = 1:length(Value.maskNames)
    outside = Value.outside{n, 1};
    for ii = 1:2
        plot(Value.potential, outside(:, ii), '.k')
    end
end
xlabel('Potential/V'); ylabel('\DeltaIntensity''');
title([expName ' \DeltaIntensity'' with Potential, Na_2SO_4'])
hold off
figPath2 = [saveRoute '\' expName '_intensityVSpotential' num2str(n) ];
saveas(img2, figPath2, 'fig')
end