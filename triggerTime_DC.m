function begin = triggerTime_DC(data)

[~, begin.pike] = max(diff(data(:, 1)));
[~, begin.CS1] = max(diff(data(:, 2)));
[~, begin.CS2] = min(diff(data(100000:end, 2)));
begin.CS2 = begin.CS2 + 100000;

begin.frame = ceil((begin.CS1 - begin.pike + 1)/10000*100) - 1;
begin.end = ceil((begin.CS2 - begin.pike + 1)/10000*100) - 1;