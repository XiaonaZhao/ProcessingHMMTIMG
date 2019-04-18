%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Frequency Response Measurement            %
%              with MATLAB Implementation              %
%                                                      %
% Author: M.Sc. Eng. Hristo Zhivomirov        11/15/14 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Amp, Ph] = freqrespmeasure(x, y)

% function: [Amp, Ph] = freqrespmeasure(x, y)
% x - first signal in the time domain
% y - second signal in the time domain
% Amp - freq. response amplitude
% Ph - freq. response phase, rad

% represent x and y as column-vectors
x = x(:);
y = y(:);

% remove the DC component
x = x - mean(x);
y = y - mean(y);

% signals length
xlen = length(x);
ylen = length(y);

% window preparation
xwin = flattopwin(xlen, 'periodic');
ywin = flattopwin(ylen, 'periodic');

% define the coherent amplification of the window
Cx = sum(xwin)/xlen;
Cy = sum(ywin)/ylen;

% fft of the first signal
X = fft(x.*xwin);

% fft of the second signal
Y = fft(y.*ywin);

% spectral peaks detection
[~, indx] = max(abs(X));
[~, indy] = max(abs(Y));

% freqeuncy response amplitude
Xamp = abs(X(indx))/xlen/Cx;
Yamp = abs(Y(indy))/ylen/Cy;
Amp = Yamp/Xamp;

% frequency response phase
Ph = angle(Y(indy)) - angle(X(indx));

end