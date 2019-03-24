tic

result  = 'G:\TaS2\20190318_TaS2_AFM_CH18SH\result\';

potential1 = (0:-0.001:-0.8)';
potential2 = (-0.799:0.001:-0.001)';
% potential1 = (0:-0.001:-1)';
% potential2 = (-0.999:0.001:-0.001)';
A2.potential = [potential1' potential2' potential1' potential2']';
clear potential1 potential2

[~, A2.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
A2.tifDir = dir(fullfile(A2.tifFile, '*.tiff'));

load('E:\TaS2\20190313_TaS2_CH18-SH\data_A2.mat')
begin = triggerTime(data, t);
A2.validDir = A2.tifDir(begin.frame:(begin.frame+length(A2.potential)));

A2.maskPath = 'E:\TaS2\20190313_TaS2_CH18-SH\MaskB3';
[~, A2.maskNames] = ReadTifFileNames(A2.maskPath);

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');
for n = 1:length(A2.maskNames)
    Mask = imread(fullfile(A2.maskPath, A2.maskNames{n}));
    mask = ~Mask;
    % span = 8;
    if sum(mask(:)) == 0
        return
    end
    
    points = ReadTifMaskPoint(A2.tifFile, A2.validDir, mask);
    
    Fs = 106;
    col = size(points, 2);
    curve = zeros(size(points));
    for ii = 1:1:col
        curve(:, ii) = lowp(points(:, ii), 1, 36, 0.1, 20, Fs); % SPR, 20;
        %         curve(:, ii) = lowp(points(:, ii), 1, 36, 0.1, 20, Fs); % A3;
        %     curve(:, ii) = lowp(points(:, ii), 44*0.01, 113*0.01, 0.1, 7, Fs); % bf, 7 ;
    end
    clear points
    
    X = (1:1:(size(curve, 1)-1))';
    dcurve = zeros(size(curve, 1)-1, size(curve, 2));
    
    for ii = 1:col
        dcurve(:, ii) = diff(curve(:, ii));
    end
    
    A2.outside{n, 1} = figSketch(dcurve);
    outside = A2.outside{n, 1};
    img = figure('color','w');
    hold on
    for ii = 1:size(outside, 2)
        plot(X, outside(:, ii), '.k')
    end
    xlabel('Frames'); ylabel('\DeltaIntensity''');
    title(['Na_2SO_4 ROI' num2str(n)])
    hold off
    figPath = [result 'A2_roi' num2str(n) ];
    saveas(img, figPath, 'fig')
    
    processBar(length(A2.maskNames), n, hwait)
    %     if length(A3.maskNames) - n <= 1
    %         waitbar(n/A3.maskNames, hwait, 'Test to be accompleted');
    %         pause(0.05);
    %     else
    %         PerStr = round(n/A3.maskNames*100);
    %         str = ['Running', num2str(PerStr), '%'];
    %         waitbar(n/A3.maskNames, hwait, str);
    %         pause(0.05);
    %     end
    
end

A2.begin = begin;

save('TaS2_0313_A2.mat', 'A2',  '-v7.3');

style = [...
    'k',... % Red 1
    'g',... % Green 2
    'b',... % Blue 3
    'c',... % Cyan 4
    'm',... % Magenta 5
    'y',... % Yellow 6
    'r',... % Black 7
    'w',... % White 8
    ];
figure('color','w');
hold on
for n = 1:length(A2.maskNames)
    outside = A2.outside{n, 1};
    for ii = 1:size(outside, 2)
        plot(A2.potential, -outside(:, ii), style(n))
    end
end
xlabel('Potential/V'); ylabel('\DeltaIntensity''');
title('All ROI in one figure, K_2SO_4')
hold off

inside = cell2mat(A2.outside);
A2.inside = zeros(size(A2.outside{1, 1}));
for n =1:length(A2.potential)
    m = n:3200:size(inside, 1);
    A2.inside(n, [1 3]) = max(inside(m, [1 3]));
    A2.inside(n, [2 4]) = min(inside(m, [2 4]));
end

toc
disp(['Running time: ',num2str(toc)]);

warndlg('Test passed', 'Warning')

MailToMe('nona1588@outlook.com');



%%
X = (1:1:4565)';
c = lowp(Y(:, 2), 10, 37, 0.1, 20, Fs);
plot(X, Y(:, 2), X, c)
xlabel('Frames'); ylabel('\DeltaIntensity''');
xlim([1550 1640])
%%
X = (1:1:(size(curve, 1)-1))';
dcurve = zeros(size(curve, 1)-1, size(curve, 2));

for n = 1:col
    dcurve(:, n) = diff(curve(:, n));
end

%%
fullfig = figure('color','w');
hold on
for n = 1:col
    %     dcurve(:, n) = diff(curve(:, n));
    plot(X, dcurve(:, n), '.k')
end
xlabel('Frames'); ylabel('\DeltaIntensity''');
hold off
% imwrite(fullfig, 'A2_thin_fullfig.tif');
% close fullfig

%%

% % this section has been made into function 'plotsketch'
% % function outside = plotSketch(dcurve)

row = size(dcurve, 1);
extreme = zeros(row ,2);
subextreme = zeros(row,2);
% temp = zeros(1, row);

for n = 1:row
    temp = dcurve(n, :);
    extreme(n, 1) = max(temp);
    extreme(n, 2) = min(temp);
    temp(temp == extreme(n, 1)) = extreme(n, 2);
    subextreme(n, 1) = max(temp);
    temp(temp == extreme(n, 2)) = subextreme(n, 1);
    subextreme(n, 2) = min(temp);
end

figure('color','w');
hold on
plot(X, extreme(:, 1), 'k', X, extreme(:, 2), 'k')
plot(X, subextreme(:, 1), 'k', X, subextreme(:, 2), 'k')
xlabel('Frames'); ylabel('\DeltaIntensity''');
xlim([1 375])
hold off

outside = [extreme subextreme];

figure('color','w');
hold on
for n = 1:size(outside, 2)
    plot(X, outside(:, n), '.k')
end
xlabel('Frames'); ylabel('\DeltaIntensity''');
xlim([1 375])
hold off

% % end

%% for A5
X = (393:1:1193)';
figure('color','w');
hold on
for n = 1:col
    plot(X, diff(curve(392:1193, n)), '.k')
end
xlim([393, 1193]);
xlabel('Frames');
ylabel('¦ÄIntensity'' (A5)');
hold off

%% for A3 A4
dCurve = zeros((size(curve, 1) - 1), size(curve, 2));
for n = 1:1:col
    dCurve(:, n) = diff(curve(:, n));
end
clear curve

for n = 1:1:col
    dCurve(:, n) = smooth(dCurve(:, n), 3);
end

row = size(dCurve, 1);
m = 1:2:row;
pdCurve = dCurve(m, :);
X = (1:1:size(pdCurve, 1))';

figure('color','w');
hold on
for n = 1:col
    plot(X, pdCurve(:, n), '.r')
end
hold off

%%
Y = [outside_A3(100:900, 3:4), outside_A4(219:1019, 3:4), outside_A5(393:1193, 3:4)];
color = [...
    'r',... % Red 1
    'g',... % Green 2
    'b',... % Blue 3
    'c',... % Cyan 4
    'm',... % Magenta 5
    'y',... % Yellow 6
    'k',... % Black 7
    'w',... % White 8
    ];
figure('color','w');
hold on
for n = 1:2:5
    plot(X, Y(:, n), color(n), X, Y(:, n+1), color(n));
end
hold off
%%
figure('color','w');
hold on
for n = 1:length(A2.maskNames)
    outside = A2.outside{n, 1};
    for ii = 1:size(outside, 2)
        plot(A2.potential(1601:3200), -outside(1601:3200, ii), '.k')
    end
end

for n = 1:length(A2.maskNames)
    outside = A2.outside{n, 1};
    for ii = 1:size(outside, 2)
        plot(A2.potential(1601:3200), -outside(1601:3200, ii), '.r')
    end
end
xlabel('Frames'); ylabel('\DeltaIntensity''');
legend('Li_2SO_4, 2nd', 'Na_2SO_4, 2nd')
hold off
%%

for n =1:3200
    m = n:3200:size(inside, 1);
    outside(n, :) = max(inside(m, :));
end
%%
figure('color','w');
hold on
for ii = 1:2
    plot(A2.potential(1601:3200), -A2.inside(1601:3200, ii), '.k')
end
for ii = 1:2
    plot(B1.potential(1601:3200), -B1.inside(1601:3200, ii), '.r')
end
xlabel('Frames'); ylabel('\DeltaIntensity''');
legend('Na_2SO_4', 'Li_2SO_4');
hold off
%%
figure('color','w');
hold on
plot(potential(2001:4000), Na_100(2001:4000), 'LineWidth', 2)
plot(potential(2001:4000), Li_100(2001:4000), 'LineWidth', 2)
plot(potential(2001:4000), Na_400(2001:4000), 'LineWidth', 2)
legend('Na_2SO_4, 100 mV/s', 'Li_2SO_4, 100 mV/s', 'Na_2SO_4, 400 mV/s')
xlabel('Potential/V'); ylabel('Current/A');
hold off