function [f, P1] = fft_P1(X, Fs)


% Fs = 800;
% X = Current(3:end, 2);
Y = fft(X);
L = length(X);
P2 = abs(Y/L);
P1 = P2(1:ceil(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:ceil(L/2))/L;
% plot(f, P1)
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')