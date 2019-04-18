% Projection aspect ratio of irregular graphics
% L = zeros(row,1);
% B = zeros(row,1);
% LvsB = zeros(row,1);
function LvsB = LengthWidth(BW)
    k = find(sum(BW));
    L = k(end) - k(1) + 1; % 
    j = find(sum(BW, 2));
    B = j(end) - j(1) + 1;
    LvsB = L/B;