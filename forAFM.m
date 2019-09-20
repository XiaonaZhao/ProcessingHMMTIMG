Z = reshape(z_data, 256, 256)';
figure('color', 'w');
imshow(Z, 'DisplayRange',[], 'InitialMagnification', 'fit');
impixelinfo

I = Z;
%%
figure('color', 'w')
subplot(221)
imshow(I, 'DisplayRange',[], 'InitialMagnification', 'fit');

%%
I2 = im2bw(mat2gray(I), 100*0.001); % for AFM, 0.1 is perfect.
subplot(222)
imshow(I2, 'DisplayRange',[], 'InitialMagnification', 'fit');
impixelinfo

%%
Me = mean(I(I2 == 0));
I0 = ones(256, 256)*Me;

%%
imLabel = bwlabel(I2);
stats = regionprops(imLabel, 'Area');
area = cat(1, stats.Area);
index = find(area == max(area));% Find the index of the smallest connected domain
I3 = ismember(imLabel, index);
% I3  = and(~I2, I3);
subplot(223)
imshow(I3, 'DisplayRange',[], 'InitialMagnification', 'fit');

%%
% imCentral = I1.*I3;
% S = sparse(imCentral);
% Ma = max(S(:));
% Mi = min(S(S>0));
% C = (S - Mi)/(Ma - Mi);
% I4 = (full(C)).*I3;
% subplot(224)
% imshow(I4, 'DisplayRange',[], 'InitialMagnification', 'fit');
% colormap jet

%%
I5 = I*(10^7);
I5(I5<0) = 0;
I5(I5>2.5) = 2.5;
figure('color', 'w')
imshow(I5, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet
colorbar

%%
I6 = I5;
I6(I5>0.35) = 0.35;
figure('color', 'w')
imshow(I6, 'DisplayRange', [], 'InitialMagnification', 'fit');
colormap jet
colorbar