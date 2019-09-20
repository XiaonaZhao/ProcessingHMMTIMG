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