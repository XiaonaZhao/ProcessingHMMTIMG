% origin
% s = zeros(row, 1);
% s0 = zeros(row, 1);
% c = zeros(row, 1);
function s = countIntensity(img, mask)

img0 = img.*mask;
% img1 = img.*mask1; % background1
% img2 = img.*mask2; % background2
% s0 = sum(img0(:))
s0 = sum(img0(:));
% bg1 = sum(img1(:));
% bg2 = sum(img2(:));
% s0 = sum(img0(:))-(bg1+bg2)/2
% c = sum(mask(:));

% I still haven't solve the problem of automatically shifting baseline.
% origin = 3.77313;

s = -s0;
% s = -(s0 - c*origin); % shifting baseline.