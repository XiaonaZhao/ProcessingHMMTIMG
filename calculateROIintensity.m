function s = calculateROIintensity(img, mask)

bitwiseMask = ~mask;
img = double(img);

% Find the position of mono- sample
centralMask = findMinDomain(mask);

% Get the outline of mono- sample
contour = edge(centralMask ,'canny'); % 'canny' would enroll some others.
figure('color', 'w'); % check figure
imshow(contour);
title('Outline');

% BWedge = bwperim(centralMask);
% figure('color', 'w'); % check figure
% imshow(BWedge);
% title('Outline');

% Normalize the intensity in the ROI
% norIMG = normalizeIntensity(img, bitwiseMask, BWedge);
norIMG = normalizeIntensity(img, bitwiseMask, contour);

% Count all normalized intensity
s = sum(norIMG(:));
