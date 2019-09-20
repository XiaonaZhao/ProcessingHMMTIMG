savepath = 'G:\TaS2\20190501_TaS2_ITO\_Result_106_std\Pics\';
%%
[row, col] = roiLocation(5);

frame = cell(67, 1);
tif0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
for ii = 1242:1308
%      tif = (double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name))) - tif0)./tif0*100;
     tif = double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name))) - tif0;
     frame{(ii-1242+1), 1} = tif((row(1):row(2)), (col(1):col(2)));
end
%%

% for ii = 17:20
    ii = 52;
    I = frame{ii, 1};
    subplot(231)
    imshow(I, 'DisplayRange',[], 'InitialMagnification', 'fit');
    
    I0 = rgb2gray(I);
    subplot(232)
    imshow(I0, 'DisplayRange',[], 'InitialMagnification', 'fit');
    
    I1 = imboxfilt(I0, 11);
    subplot(233)
    imshow(I1, 'DisplayRange',[], 'InitialMagnification', 'fit');
    
    I2 = im2bw(mat2gray(I1),259*0.001);
    subplot(234)
    imshow(I2, 'DisplayRange',[], 'InitialMagnification', 'fit');
    
    imLabel = bwlabel(~im2bw(mat2gray(imboxfilt(rgb2gray(frame{22, 1}), 13)), 242*0.001));
    stats = regionprops(imLabel, 'Area');
    area = cat(1, stats.Area);
    index = find(area == max(area));% Find the index of the smallest connected domain
    centralBW = ismember(imLabel, index);
    I3  = and(~I2, centralBW);
    subplot(235)
    imshow(I3, 'DisplayRange',[], 'InitialMagnification', 'fit');
    
    imLabel = bwlabel(~I3);% Mark each connected domain
    stats = regionprops(imLabel, 'Area');% Count the size of each connected domain
    area = cat(1, stats.Area);
    index = find(area == max(area));% Find the index of the smallest connected domain
    centralBW = ismember(imLabel, index);
    subplot(236)
    imshow(centralBW, 'DisplayRange',[], 'InitialMagnification', 'fit');
    
%     imwrite(centralBW, [savepath, '\' 'C2_' num2str(ii, '%04d'), '.tif']);
%     
%     close gcf
% end
%%
% for ii = 1:size(frame, 1)
% for ii = 21:36
    % from 21 -> 36 TaS2_20190506_C2_video_frames
    ii = 64;
    I = frame{ii, 1};
    
    subplot(221)
    imshow(I, 'DisplayRange',[], 'InitialMagnification', 'fit');
%     I0 = rgb2gray(I);
%     subplot(232)
%     imshow(I0, 'DisplayRange',[], 'InitialMagnification', 'fit');
    I1 = imboxfilt(I, 3);
    subplot(222)
    imshow(I1, 'DisplayRange',[], 'InitialMagnification', 'fit');

    I2 = im2bw(mat2gray(I1),500*0.001);
    %     I_12 = binary_otus(I_10);
    subplot(223)
    imshow(I2, 'DisplayRange',[], 'InitialMagnification', 'fit');

    imLabel = bwlabel(~I2);% Mark each connected domain
    stats = regionprops(imLabel, 'Area');% Count the size of each connected domain
    area = cat(1, stats.Area);
    index = find(area == max(area));% Find the index of the smallest connected domain
    centralBW = ismember(imLabel, index);
    subplot(224)
    imshow(~centralBW, 'DisplayRange',[], 'InitialMagnification', 'fit');
    %     I2 = im2bw(mat2gray(I1), reT);
%     imwrite(~centralBW, [savepath, '\' 'C2_' num2str(ii, '%04d'), '.tif']);
%     %     saveas(gcf,[savepath, '\' 'C2_' num2str(ii, '%04d'), '.tif']);
%     close gcf
% end
