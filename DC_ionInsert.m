% imread tiff files and cut them into PC level
tifFile = 'G:\TaS2\TaS2_20190520_ITO_AFM\D2_Na_-0-8V_25s_Pike100fps';
tifDir = dir(fullfile(tifFile, '*.tiff'));

[row, col] = ImageJroiLocation(sROI{2});
sub = cell(199, 1);
for ii = frame1:frame2
    tif = im2double(imread(fullfile(tifFile, tifDir(ii+1).name)))...
        - im2double(imread(fullfile(tifFile, tifDir(ii).name)));
    sub{ii-frame1+1} = tif((row(1):row(2)), (col(1):col(2)));
end


%% introduce the 'Mask'
Mask1 = imread('G:\TaS2\TaS2_20190520_ITO_AFM\Mask\MaskD2_roi2.tif');
a = ones(224, 1);
Mask1 = [Mask1 a];
a = ones(1, 221);
Mask1 = [Mask1; a];
nn = ones(size(Mask1));


%%

% ====== for Accumulation pic ========

% savepath = 'E:\TaS2\TaS2_20190520_ITO_AFM\Result\D2_roi2_sum2';
Viewer = zeros(199, 6);

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
    %     index = find(area == max(area));% Find the index of the smallest connected domain
    %     centralBW = ismember(imLabel, index);
    %     subplot(236)
    %     imshow(~centralBW, 'DisplayRange',[], 'InitialMagnification', 'fit');
    
    I3 = sub{ii}.*((~I2)); % !!!!!!!!!!!!!!!!!!!!!!!!****==============
    % imshow(I3, 'DisplayRange',[], 'InitialMagnification', 'fit');
    % pause(0.01);
    pic = pic+I3/199;
    
    Viewer(ii, 1) = length(area);
    if isempty(area)
        Viewer(ii, 2) = 0;
    else
        Viewer(ii, 2) = max(area);
    end
    Viewer(ii, 3) = sum(pic(:));
    
%     Marker0 = I2 | Mask1;
%     Marker = Marker + Marker0;
%     Viewer(ii, 4) = length(find(Marker < 50));
%     Viewer(ii, 5) = length(find(Marker > 20 & Marker < 120));
%     Viewer(ii, 6) = length(find(Marker >= 100));
end
figure('color','white');
plot(x, Viewer(:, 3))
 
%% Output the video
savepath = 'G:\TaS2\TaS2_20190520_ITO_AFM\Result\D2_theabove_default';
Viewer = zeros(199, 2);
video = VideoWriter('D2_theabove_default.avi'); %初始化一个avi文件
video.FrameRate = 10;
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
    %     index = find(area == max(area));% Find the index of the smallest connected domain
    %     centralBW = ismember(imLabel, index);
    %     subplot(236)
    %     imshow(~centralBW, 'DisplayRange',[], 'InitialMagnification', 'fit');
    
    I3 = nonsub.*((~I2));
    % imshow(I3, 'DisplayRange',[], 'InitialMagnification', 'fit');
    % pause(0.01);
    pic = pic+I3/199;
    imshow(pic, 'DisplayRange',[], 'InitialMagnification', 'fit');
    colormap default
    title(['t = ', num2str((ii-1)/100), ' s'], 'FontSize', 14)
    colorbar('FontSize',14)
    caxis([0 0.02])
    scalebar();
    pause(0.01);
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, [savepath, '\' 'D2_' num2str(ii, '%04d'), '.tif']);
    
    writeVideo(video, imagewd.cdata);
    
    Viewer(ii, 1) = length(area);
    if isempty(area)
        Viewer(ii, 2) = 0;
    else
        Viewer(ii, 2) = max(area);
    end
    
end

close(video);



%% in case that I forgot the HEAD

tif140 = double(imread(fullfile(tifFile, tifDir(140).name)));
tif = cell(200, 1);
for jj = 140:339
    tif0 = double(imread(fullfile(tifFile, tifDir(jj).name))) - tif140;
    tif{jj-140+1} = tif0(row(1):row(2), col(1):col(2));
end
%% debug the parameter in 'imboxfilt'
Marker = zeros(200, 3);
for ii = 1:200
    
    I = imboxfilt(tif{ii}, 13);
    I = I.*(~Mask1);
%     Marker(ii, 1) = length(find(I > 0)) + length(find(I < 0 & I > -100)) + length(find(I < -200));
    Marker(ii, 3) = 0.125*0.125*(length(find(I > 0)) + length(find(I < 0 & I > -120)));
    Marker(ii, 2) = 0.125*0.125*length(find(I <= -80 & I > -180));
    Marker(ii, 1) = 0.125*0.125*length(find(I <= -150));
end

%% ====== find the longest distance inside =====
% calculate the ion diffusion

rectXLocation = zeros(199, 4);
rectYLocation = zeros(199, 4);
RecordMaxArea = zeros(199, 1);

% pic = zeros(size(Mask1));
Marker = zeros(size(Mask1));
for ii = 1:199
    
    %         ii = 29;
    nonsub = (sub{ii}+nn)/2;
    % %%
    % BgValue = mean(nonsub(Mask1 ~= 0));
    % nonsub = (nonsub-BgValue*ones(size(nonsub)))/(1+BgValue);
    %     subplot(231)
    %     imshow((sub{ii}+nn)/2, 'DisplayRange',[], 'InitialMagnification', 'fit');
    %     pause(0.01);
    
    I = single(nonsub);
    
    [L, Centers] = imsegkmeans(I,3);
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
            RecordMaxArea(ii, 1) = 0;
            break
        end
        I2 = im2bw(mat2gray(I1), reT2*0.001);
        I2 = I2 | Mask1;
        imLabel = bwlabel(~I2);% Mark each connected domain
        stats = regionprops(imLabel, 'Area');% Count the size of each connected domain
        area = cat(1, stats.Area);
    end
    
    if isempty(area) == 0 && max(area) == 1
        RecordMaxArea(ii, 1) = 1;
    end
    
    
    if isempty(area) == 0 && max(area) > 1
        I3 = ~I2;
        RecordMaxArea(ii, 1) = max(area);
        ind = find(area == RecordMaxArea(ii, 1));%找到最大连通区域的标号
        [r, c]=find(imLabel == ind(1));
        % 获取标记物体最小外接矩形坐标点
        [rectx,recty,~,~] = minboundrect(c,r,'p');
        
        for jj = 1:4
            rectXLocation(ii, jj) = rectx(jj, 1);
            rectYLocation(ii, jj) = recty(jj, 1);
        end
    end
end

%% Matlab质心求扩散系数
ionCenter = zeros(199, 2); % the centroid coordinates
Marker = zeros(size(Mask1));
for ii = 1:199
    
    % ii = 29;
    nonsub = (sub{ii}+nn)/2;
    I = single(nonsub);
    
    [L, Centers] = imsegkmeans(I,3);
    B = labeloverlay(I,L);
    
    reT1 = 13;
    I1 = imboxfilt(I, reT1);

    reT2 = 300;
    I2 = im2bw(mat2gray(I1), reT2*0.001);
    
    I2 = I2 | Mask1;
    
    I3 = ~I2;
    imLabel = bwlabel(I3);% Mark each connected domain
    stats = regionprops(imLabel, 'Area');% Count the size of each connected domain
    area = cat(1, stats.Area);
    while isempty(area)
        reT2 = reT2+50;
        if reT2 > 450
            ionCenter(ii, 1) = ionCenter(ii-1, 1); 
            ionCenter(ii, 2) = ionCenter(ii-1, 2);
            break
        end
        I2 = im2bw(mat2gray(I1), reT2*0.001);
        I2 = I2 | Mask1;
        I3 = ~I2;
        imLabel = bwlabel(I3);% Mark each connected domain
        stats = regionprops(imLabel, 'Area');% Count the size of each connected domain
        area = cat(1, stats.Area);
    end
    
    if isempty(area) == 0
        ind = find(area == max(area));
        I3(find(imLabel ~= ind(1))) = 0;
        I4 = nonsub.*(I3);
        [ionCenter(ii, 1), ionCenter(ii, 2)] = oCenter(I4); % get the centroid of the maximun connected domain
    end
    
end

%% calculate the D
D = zeros(198, 1);
t = 0.01;
m = 0;
for ii =1:198
    pixel = sqrt((ionCenter(ii+1, 1)-ionCenter(ii, 1))^2+(ionCenter(ii+1, 2)-ionCenter(ii, 2))^2);
    distance = pixel*125*(10^(-9));
    D(ii, 1) = 4*(distance^2)/(pi^2*t);
    
    if distance == 0
        m = m + 1;
        t = t + t*m;
    else
        t = 0.01;
        m = 0;
    end

end
%% Matlab质心求扩散系数 average 5 span
span = 10;
ionCenter = zeros(fix(200/span), 2); % the centroid coordinates
Marker = zeros(size(Mask1));
for ii = 1:span:199
    
    % ii = 29;
    nonsub = zeros(size(nn));
    for jj = 0:span-1
        nonsub = nonsub + (sub{ii+jj}+nn)/2;
    end
    
%     nonsub = nonsub/span;
    I = single(nonsub);
    
    [L, Centers] = imsegkmeans(I, 3);
    B = labeloverlay(I,L);
    
    reT1 = 13;
    I1 = imboxfilt(I, reT1);

    reT2 = 300;
    I2 = im2bw(mat2gray(I1), reT2*0.001);
    
    I2 = I2 | Mask1;
    
    I3 = ~I2;
    imLabel = bwlabel(I3);% Mark each connected domain
    stats = regionprops(imLabel, 'Area');% Count the size of each connected domain
    area = cat(1, stats.Area);
    while isempty(area)
        reT2 = reT2+50;
        if reT2 > 450
            if ii == 1
                ionCenter(ceil(ii/span), 1) = 52;
                ionCenter(ceil(ii/span), 2) = 141;
            else
                ionCenter(ceil(ii/span), 1) = ionCenter(ceil(ii/span)-1, 1);
                ionCenter(ceil(ii/span), 2) = ionCenter(ceil(ii/span)-1, 2);
            end
            break
        end
        I2 = im2bw(mat2gray(I1), reT2*0.001);
        I2 = I2 | Mask1;
        I3 = ~I2;
        imLabel = bwlabel(I3);% Mark each connected domain
        stats = regionprops(imLabel, 'Area');% Count the size of each connected domain
        area = cat(1, stats.Area);
    end
    
    if isempty(area) == 0
        ind = find(area == max(area));
        I3(find(imLabel ~= ind(1))) = 0;
        I4 = nonsub.*(I3);
        [ionCenter(ceil(ii/span), 1), ionCenter(ceil(ii/span), 2)] = oCenter(I4); % get the centroid of the maximun connected domain
    end
    
end

D = zeros(19, 1);
t = 0.01;
m = 0;
for ii =1:fix(200/span)-1
    pixel = sqrt((ionCenter(ii+1, 1)-ionCenter(ii, 1))^2+(ionCenter(ii+1, 2)-ionCenter(ii, 2))^2);
    distance = pixel*125*(10^(-9));
    D(ii, 1) = 4*(distance^2)/(pi^2*t);
    
    if distance == 0
        m = m + 1;
        t = t + t*m;
    else
        t = 0.01;
        m = 0;
    end

end
time = 0.01*span*(1:length(D))';
% figure('color', 'w');
plot(time, D)