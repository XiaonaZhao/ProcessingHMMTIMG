function outside = figSketch(dcurve)

row = size(dcurve, 1);
extreme = zeros(row ,2);
subextreme = zeros(row,2);
% temp = zeros(1, row);

for n = 1:row
    temp = dcurve(n, :);
    extreme(n, 1) = max(temp);
    extreme(n, 2) = min(temp);
    temp(temp == extreme(n, 1)) = extreme(n, 2);
    subextreme(n, 1) = max(temp);
    temp(temp == extreme(n, 2)) = subextreme(n, 1);
    subextreme(n, 2) = min(temp);
end

outside = [extreme subextreme];

end