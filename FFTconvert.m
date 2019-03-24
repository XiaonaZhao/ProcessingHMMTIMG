function output = FFTconvert(img, mask)

img1 = fft2(img);
img2 = fftshift(img1);
img2 = img2.*mask;
img1 = ifftshift(img2);
output = ifft2(img1);