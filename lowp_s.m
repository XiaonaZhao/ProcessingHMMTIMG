function output = lowp_s(input, Fs)

frame = size(input, 3);
output = zeros(size(input));
for ii = 1:size(input, 1)
    parfor jj = 1:size(input, 2)
        temp = input(ii, jj, :);
        temp = reshape(temp, frame, 1);
        temp = lowp(temp, 4, 16, 33*0.01, 20, Fs); % Set Fs = 5 Hz;
                temp = highpass(temp, 4.5, Fs);
%         temp = lowp(temp, 9, 18, 33*0.01, 20, Fs); % Set Fs = 10 Hz;
%         temp = highpass(temp, 9.5, Fs);
        output(ii, jj, :) = temp;
    end
end
end
