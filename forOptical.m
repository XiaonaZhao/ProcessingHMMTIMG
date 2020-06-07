% begin = triggerTime_DC(data);
begin = triggerTime_AC(data, t, 100);

%%
[~, Value.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));
I0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
% I = double(imread(fullfile(Value.tifFile, Value.tifDir(begin.end-30).name)));
I = double(imread(fullfile(Value.tifFile, Value.tifDir(begin.peak2).name)));
I1 = I./I0;
figure('color', 'w')
imshow(I1, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap(flipud(jet))
colorbar
impixelinfo

%%
figure('color', 'w')
subplot(221)
imshow(I, 'DisplayRange', [], 'InitialMagnification', 'fit');


I2 = ~im2bw(mat2gray(I), 291*0.001); % for optical
subplot(222)
imshow(I2, 'DisplayRange', [], 'InitialMagnification', 'fit');


imLabel = bwlabel(I2);
stats = regionprops(imLabel, 'Area');
area = cat(1, stats.Area);
index = find(area == max(area));% Find the index of the smallest connected domain
I3 = ismember(imLabel, index);
% I3  = and(~I2, I3);
subplot(223)
imshow(I3, 'DisplayRange', [], 'InitialMagnification', 'fit');


imCentral = I1.*I3;
S = sparse(imCentral);
Ma = max(S(:));
Mi = min(S(S>0));
C = (S - Mi)/(Ma - Mi);
I4 = (full(C)).*I3; 
subplot(224)
imshow(I4, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet

figure('color', 'w')
imshow(I4, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet
colorbar

%%
I5 = (I-I0)./I0.*I3;
% I5(I5<0) = 0;
figure('color', 'w')
imshow(I5, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet
colorbar


I6 = I5;
I6(I5>0.35) = 0.35;
figure('color', 'w')
imshow(I6, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet
colorbar
%%
[cstrFilenames, cstrPathname] = uigetfile(...
    {'*.*',  'All Files (*.*)';...
    '*.zip',  'Zip-files (*.zip)';...
    '*.roi',  'ROI (*.roi)'...
    },'Pick a .roi imageJ file');
[sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
ROI = sROI{1, 1};
[row, col] = ImageJroiLocation(ROI);

%%
L = 1024;
prefix = ('G:\MoS2\MoS2_0919_0802\_Result\A_sequance\');
d = sortObj(dir([prefix, '*.mat']));
intensity = zeros(L, 6);
% for jj = 1:6
jj = 4;
    load([prefix, d(jj).name]);
    mask = ~imread(expTab(jj).roiMask);
    tif0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
    exp = cell(L, 1);
    
    for ii = 1:L
        temp = double(imread(fullfile(Value.tifFile, Value.validDir(ii).name)));
        exp{ii, 1} = (temp - tif0)./tif0;
%         intensity(ii,jj) = ROImean(temp, mask);
        %     tif{ii, 1} = temp(row(1):row(2), col(1):col(2));
    end
% end
%%
Fs = 100; % normal hamamatsu camera
for ii = [1 2 4 5 6]
    Y = intensity(:, ii);
    [f1(:, ii), P1(:, ii)] = fft_P1(Y, Fs);
%     figure('color','w');
    plot(f1(:, ii), 1000*P1(:, ii));
    hold on
    xlim([0, 50]);
    ylim([0 8]);
    xlabel('f (Hz)','fontsize',12)
    ylabel('|P(f)|','fontsize',12)
%     l = num2str(ii);
%     legend(l)
    % hold off
end
hold off