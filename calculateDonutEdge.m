function [insideEdge, outsideEdge] = calculateDonutEdge(img, mask)
% img is double one

[insideRing, outsideRing] = detectDonutEdge(mask); 
insideEdge = ROImean(img, insideRing);
outsideEdge = ROImean(img, outsideRing);
end

% BWedge = bwperim(centralMask);
% figure('color', 'w'); % check figure
% imshow(BWedge);
% title('Outline');

% Normalize the intensity in the ROI
% norIMG = normalizeIntensity(img, bitwiseMask, BWedge);
% norIMG = normalizeIntensity(img, bitwiseMask, contour);

% Count all normalized intensity
% s = sum(norIMG(:));