function begin = triggerTime_PEIM(data, t, Fs)

[~, begin.pike] = max(diff(data(:, 1)));

plot(t, data(:, 2))
xlim([4 5])

[x, ~] = ginput(1);
begin.CS1 = x*10000;

begin.frame = ceil((begin.CS1 - begin.pike + 1)/10000*Fs);

end