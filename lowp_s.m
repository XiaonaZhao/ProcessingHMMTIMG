function output = lowp_s(input, Fs)

frame = size(input, 3);
output = zeros(size(input));
for ii = 1:size(input, 1)
    for jj = 1:size(input, 2)
        temp = input(ii, jj, :);
        temp = reshape(temp, frame, 1);
        temp = lowp(temp, 4, 16, 33*0.01, 20, Fs);
        temp = highpass(temp, 4.5, Fs);
        output(ii, jj, :) = temp;
    end
end
