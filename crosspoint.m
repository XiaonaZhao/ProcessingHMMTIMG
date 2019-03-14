function [px, py] = crosspoint(x1, y1, x2, y2)
% y = kx + b
p1 = polyfit(x1, y1, 1);
p2 = polyfit(x2, y2, 1);
k1 = p1(1); b1 = p1(2);
k2 = p2(1); b2 = p2(2);

if (k1 == k2) && (b1 == b2)
    disp('The lines are on the same line.');
elseif (k1 == k2) && (b1 ~= b2)
    disp('The is no crossed point.');
else
    px = (b2 - b1)./(k1 - k2);
    py = k1*px + b1;
end
end