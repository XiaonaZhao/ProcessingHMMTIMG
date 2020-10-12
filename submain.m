% Raw data is from the sampling area subtracting the background.

% 1 raw line
% 2 lowp
% 3 windowSize % if it is necessary
% 4 drift baseline
% ues detrend
% 5 minima and its location to value1
% baseline value on the minima location to value2
% 6 minima substracts baseline (output = value2 - value1)
% 7 repeat step 1 to 6 and list output into a column vector
% actually I looped the steps 1 to 6
% 8 plot the value to col

% NOTICE:
% The submain is used to find out the diffusion boundary with the method
% from Hypo. The method is based on the attenuator of signal-to-noise ratio
% (SNR).

% Created on 16 January 2019



%% -- raw line

% from dataINxlsx.m
X = sheet(:, 1);
Y = sheet(:, (2:end));

col = size(Y, 2); % 2 returns the num of col
Yfit = zeros(size(Y));
for n = 1:col
    %     Yfit(:, n) = log(abs(Y(:, n)));
    Yfit(:, n) = dataFit(Y(:, n));
end

figure('color','w');
for n = 1:col
    plot(X, Y(:, n)); % get raw lines
    %     plot(X, Yfit(:, n)); % get fitted lines
    hold on
end

%% -- have a view of signal frequency spectrum

% - signal frequency spectrum
% Fs = 50; % low-speed hamamatsu camera
Fs = 100; % normal hamamatsu camera
% Fs = 200; % normal hamamatsu camera with darker field
% Fs = 1600; % high-speed hamamatsu camera
for n = 1:col
    [f1, P1] = fft_P1(Y(:, n), Fs);
    figure('color','w');
    plot(f1, P1);
    % hold on
    xlim([1, 15]);
    xlabel('f (Hz)','fontsize',10)
    ylabel('|P(f)|','fontsize',10)
    l = num2str(n);
    legend(l)
    % hold off
end

%% -- lowp

% Chebyshev low pass filter
filterY = zeros(size(Y));
for n = 1:col
    %     filterY(:, n) = lowp(Y(:, n), 5, 29, 0.1, 20, Fs); % Fs = 1600
    filterY(:, n) = lowp(Y(:, n), 4, 30, 0.1, 20, Fs); % Fs = 100
end

%% Check the lowp filter
% i = 2;
figure('color','w');
% plot(X, Y(:, i))
plot(X(4:end), b(3:end))
% plot(number, Y)
hold on
% filterY(:, i) = lowp(Y(:, i), 5, 29, 0.1, 20, Fs);% Fs = 1600
% filterY = lowp(Y, 5*0.1, 92*0.1, 7*0.01, 9, Fs); % Fs = 50
% filterY = lowp(Y, 4, 21, 0.1, 20, Fs); % Fs = 200
filterB = lowp(b, 2, 21, 32*0.1, 188*0.1, Fs); % Fs = 200
% plot(X, filterY(:, i))
plot(X(4:end), filterB(3:end))
% plot(number, filterY)
% xlim([0, 25600]) % Fs = 1600
% xlim([1, 800])
hold off
%% DeltaI2I_multiVSmonolayer_1v1_1
% CurrentDensity_multiVSmonolayer_1v1_1

% figure('color', 'w');
X = (0:400)';
ii = 404;
multi_i = lowp(multi/10, 3*0.1, 15, 0.01, 20, 100);
single_i = smooth(single/10, 10);
plot(X, multi_i(ii:ii+400)-multi_i(ii), X, single_i(ii:ii+400)-single_i(ii))
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
set(gca, 'linewidth', 1, 'FontSize', 16)
xlabel('Potential (V vs. Ag/AgCl)'); 
ylabel('Current density (a.u.)');

multi_c = intensity2current(multi_i, 850);
single_c = intensity2current(single_i, 850);
%%
t1 = 404;
t2 = 404;
plot(potential(1:400), -multi_c(t1:t1+399), potential(1:400), -single_c(t2:t2+399))
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
set(gca, 'linewidth', 1, 'FontSize', 16)
xlabel('Potential (V vs. Ag/AgCl)'); 
ylabel('Current density (a.u.)');

%% re-check
% Fs = 1600; % high-speed hamamatsu camera
% [f1, P1] = fft_P1(filterY(:, i), Fs);
% figure('color','w');
% plot(f1, P1);
% % hold on
% xlim([1, 20]);
% % ylim([0, 6]);
% xlabel('f (Hz)','fontsize',10)
% ylabel('|P(f)|','fontsize',10)
% l = num2str(i);
% legend(l)
% %% re-filter
% i = 1;
% % figure('color','w');
% plot(X, filterY(:, i))
% hold on
% filter1Y(:, i) = lowp(filterY(:, i), 10*0.1, 18*0.1, 10*0.1, 20*0.1, Fs);
% plot(X, filter1Y(:, i))
% xlim([0, 25600])
% hold off

%% windowSize filter

windowSize = 4.5*100;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
winfilterY = zeros(size(Y));
for n =1:col
    winfilterY(:, n) = filter(b, a, filterY(:, n));
    figure('color','w');
    plot(X, filterY(:, n))
    hold on
    plot(X, winfilterY(:, n))
    l = num2str(n);
    legend('Input Data','Filtered Data')
    title(l)
    hold off
end

% %% Check the windowSize filter
% windowSize = 9*10;
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% i = 1;
% % figure('color','w');
% plot(X, filterY(:, i))
% hold on
% winfilterY(:, 1) = filter(b, a, filterY(:, 1));
% plot(X, winfilterY(:, i))
% xlim([0, 25600])
% hold off

%% -- drift baseline

% save comparing figures
reName = 'detrend_';
INpath = 'E:\20190116_MoS2_CH18-SH_0108\detrend figures';

% x = [0 13182 25601]';
% y = [-0.0003 -47.19 -64.74]';
detrendY = zeros(size(Y));
for n = 1:col
    figure('color','w')
    plot(X, filterY(:, n))
    hold on
    [x,y] = ginput;
    % ginput gathers an unlimited number of points until you press the RETURN
    % key.
    p = polyfit(x,y,1);
    f = polyval(p,X);
    % figure('color','w');
    detrendY(:, n) = filterY(:, n) - f;
    % figure('color','w');
    plot(X, detrendY(:, n));
    legend('data','linear fit')
    hold off
    
    sName = num2str(n);
    INfilename = [reName sName '.fig'];
    SavePath = fullfile(INpath, INfilename);
    saveas(gcf, SavePath);
    close(gcf)
end

%% -- minima and its location to value1
% baseline value on the minima location to value2

minimaY = zeros(col,1);
value1 = zeros(col,1);
value2 = zeros(col,1);
value = zeros(col,1);
for n = 1:col
    [minimaY(n), loc] = min(detrendY(:, n));
    value1(n) = detrendY(loc, n);
    value2(n) = f(loc); % local baseline value
    value(n) = value2(n) - value1(n); % distance between signal and noise
end
%% 20200425
L = length(Value.tifName);
% sROI = Value.sROI;
tifPath = Value.tifPath;
tifName = Value.tifName;
tif1 = im2double(imread(fullfile(tifPath, tifName{1})));

average = zeros(L, length(sROI));
for jj = 1:L
    tif = im2double(imread(fullfile(tifPath, tifName{jj}))) - tif1;
    
    for ii = 1:length(sROI)
        [row, col] = ImageJroiLocation(sROI{ii});
        tif0 = tif(row(1):row(2), col(1):col(2));
        average(jj, ii) = mean(tif0(:));
    end
    
end
%%
for jj = 140:339
    tif0 = im2double(imread(fullfile(tifFile, tifDir(jj).name))) - tif140;
    tif{jj-140+1} = tif0(row(1):row(2), col(1):col(2));
end
%%
for jj = 1:199
    sub{jj} = tif{jj+1} - tif{jj};
end

%%
savepath = 'E:\TaS2\TaS2_20190520_ITO_AFM\Result\D2_roi2_matlabPic';
Viewer = zeros(199, 2);
video = VideoWriter('D2_roi2_sum2.avi'); %初始化一个avi文件
video.FrameRate = 24;
open(video);
pic = zeros(size(Mask1));
for ii = 1:199

%     ii = 44;
    nonsub = (sub{ii}+nn)/2;
    % %%
    % BgValue = mean(nonsub(Mask1 ~= 0));
    % nonsub = (nonsub-BgValue*ones(size(nonsub)))/(1+BgValue);
%     subplot(231)
%     imshow((sub{ii}+nn)/2, 'DisplayRange',[], 'InitialMagnification', 'fit');
%     pause(0.01);
    
    I = single(nonsub);
    
    [L,Centers] = imsegkmeans(I,3);
    B = labeloverlay(I,L);
%     subplot(232)
%     imshow(B, 'DisplayRange',[], 'InitialMagnification', 'fit')
%     title(['Labeled Image ', num2str(ii)])
%     pause(0.01);

    reT1 = 13;
    I1 = imboxfilt(I, reT1);
%     subplot(233)
%     imshow(I1, 'DisplayRange',[], 'InitialMagnification', 'fit');
%     pause(0.01);
    % I = uint8(sub{ii}); % No distinction at all
    % [~, reT] = binary_otus(I1);
    reT2 = 300;
    I2 = im2bw(mat2gray(I1), reT2*0.001);
    % I2 = binary_iterate(I1);
    % I2 = binary_bernsen(I1);
%     subplot(234)
%     imshow(I2, 'DisplayRange',[], 'InitialMagnification', 'fit');
%     pause(0.01);
    
    I2 = I2 | Mask1;
%     subplot(235)
%     imshow(I2, 'DisplayRange',[], 'InitialMagnification', 'fit');
%     pause(0.01);
    
    imLabel = bwlabel(~I2);% Mark each connected domain
    stats = regionprops(imLabel, 'Area');% Count the size of each connected domain
    area = cat(1, stats.Area);
    while isempty(area)
        reT2 = reT2+50;
        if reT2 > 450
            break
        end
        
        I2 = im2bw(mat2gray(I1), reT2*0.001);
%         subplot(234)
%         imshow(I2, 'DisplayRange',[], 'InitialMagnification', 'fit');
        I2 = I2 | Mask1;
%         subplot(235)
%         imshow(I2, 'DisplayRange',[], 'InitialMagnification', 'fit');
        imLabel = bwlabel(~I2);% Mark each connected domain
        stats = regionprops(imLabel, 'Area');% Count the size of each connected domain
        area = cat(1, stats.Area);
    end
    index = find(area == max(area));% Find the index of the smallest connected domain
    centralBW = ismember(imLabel, index);
    subplot(236)
    imshow(~centralBW, 'DisplayRange',[], 'InitialMagnification', 'fit');

I3 = nonsub.*((~I2));
imshow(I3, 'DisplayRange',[], 'InitialMagnification', 'fit');
% pause(0.01);
% pic = pic+I3/199;
% imshow(pic, 'DisplayRange',[], 'InitialMagnification', 'fit');
% % imshow(pic, 'DisplayRange',[], 'InitialMagnification', 'fit');
% title(['t = ', num2str((ii-1)/100), ' s'])
% colorbar
% caxis([0 0.06])
% scalebar();
pause(0.01);
% imagewd = getframe(gcf);
% imwrite(imagewd.cdata, [savepath, '\' 'D2_' num2str(ii, '%04d'), '.tif']);
imwrite(I3, [savepath, '\' 'D2_' num2str(ii, '%04d'), '.tif']);

writeVideo(video, imagewd.cdata);    

    Viewer(ii, 1) = length(area);
    if isempty(area)
        Viewer(ii, 2) = 0;
    else
        Viewer(ii, 2) = max(area);
    end
    
end

close(video);

%%
figure('color','white');
for ii = 22:23
    plot(X, average(:, ii))
    hold on
end


%% -- plot and multiple linear regression
N = (1:col)';
figure('color','w');
plot(N, value, 'o');

% IPOINTS = findchangepts(..., 'Statistic', STAT) specifies the type of
%   change to detect.  The default is 'mean':
%      'mean'   detect changes in mean
%      'rms'    detect changes in root-mean-square level
%      'std'    detect changes in standard deviation
%      'linear' detect changes in mean and slope
iN = findchangepts(value,'Statistic','linear','MinThreshold',var(value));
% iN = findchangepts(Y,'MaxNumChanges',1,'Statistic','linear')


% joinpointValue = lsq_lut_piecewise(N, value, 2);
% figure('color','w');
% plot(N, value, '.', N, joinpointValue, '+-')
% legend('experimental data (N, value)','LUT points (N, joinpointValue)')
% title('Piecewise 1-D look-up table least square estimation')
%% plot joinpoint and power fitting
% use the classification of iN to ployfit subsections, and plot all
% elements in one figure.

a1 = N(1:3); b1 = value(1:3); c1 = polyfit(a1,b1,1); d1 = polyval(c1, a1);
a2 = N(4:10); b2 = value(4:10);
c2 = polyfit(a2,b2,1); d2 = polyval(c2, a2);
a3 = fit(N, value, 'power1'); % cfit
figure('color', 'w');
plot(a3, N, value, 'o');
hold on
plot(a1, d1, 'k'); plot(a2, d2, 'k');
% plot(a3, N, value,'r');
hold off
%%

% h=plot(a,b,'k-o','Markersize',7,'Markerface','white','linewidth',1.0);
xlabel('Modulation frequency (Hz)','fontsize',12)
ylabel('Amplitude (\DeltaT/T * 10^-^3)','fontsize',12)
% Legend('a','b',0)
% hh = findobj('tag','legend');   %|
% set(hh,'fontsize',10)         %| 设置legend字号大小
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);        %| 设置图形线宽
set(gca, 'linewidth', 1.5)      %| 设置图形外边框的线宽1.5
% set(gca,'box','off')          %| 去图形外筐
% %| 设置坐标轴字号12 ，斜体，正
% set(gca,'fontsize',12,'fontweight','normal','fontangle','italic')
% %| 设置x轴labal字体为斜体,黑体，字号12
% set(get(gca,'xlabel'),'fontangle','italic','fontweight','bold', 'fontsize',12)
% %| 设置y轴labal字体为斜体，非黑体，字号12
% set(get(gca,'ylabel'),'fontangle','italic','fontweight','normal', 'fontsize',12)