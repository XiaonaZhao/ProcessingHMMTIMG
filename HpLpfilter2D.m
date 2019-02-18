function g = HpLpfilter2D(f, hpth, lpth)

PQ = paddedsize(size(f)); % high pass
D0 = hpth*PQ(1);
HBW = hpfilter('btw', PQ(1), PQ(2), D0, 2);
H = 0.5 + 2*HBW; % (0.5, 2.0)
ghf = dftfilt(f, H);
sharp = ghf + f;

PQ = paddedsize(size(sharp)); % lowpass
[U, V] = dftuv(PQ(1), PQ(2));
D0 = lpth*PQ(2);
% F = fft2(tif1, PQ(1), PQ(2));
H = exp(-(U.^2 + V.^2)/(2*(D0^2)));
g = dftfilt(sharp, H);

% figure, imshow(fftshift(H), [])
% figure, imshow(log(1 + abs(fftshift(F))), [])

figure
subplot(1,2,1), imshow(f), title('original image');
subplot(1,2,2),imshow(g, []), title('HpLpfiltered');
impixelinfo