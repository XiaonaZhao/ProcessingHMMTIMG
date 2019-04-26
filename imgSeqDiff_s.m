function [value, timer] = imgSeqDiff_s(input)

frame = size(input, 3);
row = size(input, 1);
col = size(input, 2);
value = zeros(row, col);
timer = zeros(row, col);
for ii = 1:row
    parfor jj = 1:col
        temp = input(ii, jj, :);
        temp = reshape(temp, frame, 1);
        [value(ii, jj), timer(ii, jj)] = max(abs(temp));
    end
end
end