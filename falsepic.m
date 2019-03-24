% input images with each 10 skips
for n = 1:(floor(row/10)+1)
    m = 10*(n-1)+1;
    tif{n, 1} = n;
    tif{n, 2} = double(imread(fullfile(tifFolder, tifNames{m}))) - tif0;
end

%
for n = 1:row
    tif{n, 1} = n;
    tif{n, 2} = double(imread(fullfile(tifFolder, tifNames{n}))) - tif0;
end
% Create warning dialog box
warndlg('Test passed','Warning')

%
f =  tif{67, 2};
F = fft2(f);
S = fftshift(log(1+abs(F)));
S = gscale(S, 'full16');

h = fspecial('sobel');
PQ = paddedsize(size(F));
H = freqz2(h,PQ(1),PQ(2));
H1 = ifftshift(H);

gs = imfilter(double(f),h);
gf = dftfilt(f,H1);

figure('color', 'w');
subplot(3,2,1);imshow(f);title('»Ò¶È¼¶Í¼Ïñ');
subplot(3,2,2);imshow(S);title('¸µÀïÒ¶ÆµÆ×');
subplot(3,2,3);imshow(abs(H),[ ]);title('ÆµÓòÂË²¨Æ÷:´¹Ö±sobleÑÚÄ£');
subplot(3,2,4);imshow(abs(H1),[ ]);title('ÆµÓò¾­fftshift´¦ÀíºóµÄÂË²¨Æ÷');
subplot(3,2,5);imshow(gs,[ ]);title('¿Õ¼äÓòÂË²¨Æ÷:´¹Ö±sobleÑÚÄ£');
subplot(3,2,6);imshow(gf,[ ]);title('¿Õ¼äÓò¾­fftshift´¦ÀíºóµÄÂË²¨Æ÷');

% lowpass filter
PQ = paddedsize(size(tif1));
[U, V] = dftuv(PQ(1), PQ(2));
Th = 0.02;
D0 = Th*PQ(2);
F = fft2(tif1, PQ(1), PQ(2));
H = exp(-(U.^2 + V.^2)/(2*(D0^2)));
g = dftfilt(tif1, H);

figure, imshow(fftshift(H), [])

figure, imshow(log(1 + abs(fftshift(tif1))), [])

figure, imshow(g, [])

% highpass filter 'gaussian'
PQ = paddedsize(size(tif1));
Th = 0.01;
D0 = Th*PQ(1);
H = hpfilter('btw', PQ(1), PQ(2), D0);
g = dftfilt(tif1, H);

figure, imshow(g, [])

% highpass filter 'btw'
PQ = paddedsize(size(f));
Th = 0.01;
D0 = Th*PQ(1);
HBW = hpfilter('btw', PQ(1), PQ(2), D0, 2);
H = 0.5 + 2*HBW;
gbw = dftfilt(f, HBW);
% gbw = gscale(gbw, 'full16');
figure, imshow(gbw, [])

ghf = dftfilt(f, H);
% ghf = gscale(ghf, 'full16');
figure, imshow(ghf, [])

sharp = ghf + f;
figure, imshow(sharp, [])

ghe = histeq(ghf, 256);
figure, imshow(ghe, [])

%
impixelinfo

% fft2
f_1 =  tif{67, 2};
%%
F_1 = fft2(f_1);
S_1 = fftshift(log(1+abs(F_1)));
S_1 = gscale(S_1, 'full16');

h_1 = fspecial('sobel');
PQ_1 = paddedsize(size(F_1));
H_1 = freqz2(h_1, PQ_1(1),PQ_1(2));
H1_1 = ifftshift(H_1);

gs_1 = imfilter(double(f_1),h_1);
gf_1 = dftfilt(f_1,H1_1);

figure('color', 'w');
subplot(3,2,1);imshow(f_1);title('»Ò¶È¼¶Í¼Ïñ_1');
subplot(3,2,2);imshow(S_1);title('¸µÀïÒ¶ÆµÆ×_1');
subplot(3,2,3);imshow(abs(H_1),[ ]);title('ÆµÓòÂË²¨Æ÷:´¹Ö±sobleÑÚÄ£_1');
subplot(3,2,4);imshow(abs(H1_1),[ ]);title('ÆµÓò¾­fftshift´¦ÀíºóµÄÂË²¨Æ÷_1');
subplot(3,2,5);imshow(gs_1,[ ]);title('¿Õ¼äÓòÂË²¨Æ÷:´¹Ö±sobleÑÚÄ£_1');
subplot(3,2,6);imshow(gf_1,[ ]);title('¿Õ¼äÓò¾­fftshift´¦ÀíºóµÄÂË²¨Æ÷_1');
%%

% hist
figure, imhist(uint16(f_1)), ylim('auto');
g_1 = histeq(uint16(f_1), 65536);
figure
subplot(1,2,1); imshow(g_1);
subplot(1,2,2); imhist(g_1), ylim('auto');

hnorm = imhist(uint16(f_1))./numel(uint16(f_1));
cdf = cumsum(hnorm);

x = linspace(0, 1, 1024);
plot(x, cdf);


figure('color', 'w');
subplot(2,2,1); imshow(f_1);
subplot(2,2,2); imshow();

%% edge detected
[g_canny_default, tc] = edge(g, 'canny');
g_canny_best = edge(g, 'canny', [30*0.001 30*0.010], 18);

figure
subplot(2,2,1)
imshow(f), title('original image');
subplot(2,2,2)
imshow(g, []), title('HpLpfiltered');
subplot(2,2,3)
imshow(K), title('Adaptive Filtered');
subplot(2,2,4)
imshow(g_canny_best) % just fine
% xlim([492 597])
% ylim([541 647])
%%

[g_sobel_default, ts] = edge(g, 'sobel');
g_sobel_best = edge(g, 'sobel', 5.5); % ts = 5.6237

% figure
imshow(g_sobel_best) % poor connectivity

%%
% Adaptive Filter
K = mywiener2(g, [5 5]);

figure
subplot(1,3,1), imshow(f), title('original image');
% subplot(1,3,2), imshow(f), title('noised image');
subplot(1,3,2), imshow(g, []), title('HpLpfiltered');
subplot(1,3,3), imshow(K), title('Adaptive Filtered');

impixelinfo

%% LMS stochastic gradient descent

%% 
fields = {'number', 'darkest', 'reduced', 'disc', 'tifFile'};
tifFile = strcat(maskTab{1,2},maskTab{1,3});

%%
expStruct(7).disc = "E:\";
tifFile = convertStringsToChars(strcat(expStruct(7).disc, expStruct(7).index));
[tifFolder, tifNames] = ReadTifFileNames(tifFile);

tif0 = double(imread(fullfile(tifFolder, tifNames{1})));
tif515 = double(imread(fullfile(tifFolder, tifNames{513}))) - tif0;

redCon515 = (-tif515)*10/(-filterCurve(401));
localSums = imboxfilt(redCon515, 11);
localSums0 = localSums;
localSums = abs(localSums);
localSums(localSums > 10) = 10;
figure('color', 'w');
imshow(localSums, 'DisplayRange',[]);
colormap jet 
map=colormap('jet');
colorbar;
h=colorbar;
set(get(h,'title'),'string','mM');
