% tic
% 
% result  = 'E:\20190306_MoS2_1123_CH18SH\result\';


%% 
varMat = load('F:\MoS2_final\MoS2_0919_0802\_Timer\A7_data.mat');
Fs = 100;
begin = triggerTime_MoS2(varMat.data, varMat.t, Fs);
expTab(6).begin = begin;

%%
load('G:\TaS2\TaS2_0909_S_Aux60\_Result\expTab.mat')

hwait = waitbar(0, 'Please wait for the test >>>>>>>>');
for m = 1:size(expTab, 1)
    expName = expTab(m).expName;
    tifPath = expTab(m).tifPath;
    mask = expTab(m).roiMask;
    begin = expTab(m).begin; 
    saveRoute = expTab(m).saveRoute;
%     Fs = 100;
    
     A = MoS2_charging(expName, tifPath, mask, begin, saveRoute);
    
    disp(['The average amplitude of ' expName ' is about ' A '.']);
    processBar(size(expTab, 1), m, hwait)
    
end
delete(hwait);




%% import image Sequence
% L = 1024;
% 
% % [~, Value.tifFile] = uigetfile('*.tiff', 'Multiselect', 'on', 'Read tif Folder');
% % Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));
% % Value.tif0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
% 
% Value.tifFile = tifPath;
% Value.tifDir = dir(fullfile(Value.tifFile, '*.tiff'));
% Value.validDir = Value.tifDir(begin.frame:(begin.frame+length(Value.potential)));
% Value.begin = begin;
% 
% tif0 = double(imread(fullfile(Value.tifFile, Value.tifDir(1).name)));
% for ii = Value.begin.frame:(Value.begin.frame+L)
%     tif  = (double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name))) - tif0)./tif0;
% end
% 
% 
% %% 
% tif_A_peak = zeros(size(tif0, 1), size(tif0, 2));
% 
% tif_3D = zeros(480, 640, 1024);
% for ii = 1:L
%     tif_3D(:, :, ii) = Value.tif{ii, 1};
% end
% 
% for ii = 1:size(tif0, 1)
%     for jj = 1:size(tif0, 2)
%         X = reshape(tif_3D(ii, jj, :), L, 1);
%         Y = fft(X);
%         P2 = abs(Y/L);
%         tif_A_peak(ii, jj) = max(2*P2(2:end-1));
%     end
% end
% figure
% imshow(tif_A_peak, 'DisplayRange', [], 'InitialMagnification', 'fit');
% colormap jet
% 
% 
% 
% 
% %% Total partly sampling analysis vs. bg
% Value.roiMask = 'I:\20190306_MoS2_1123_CH18SH\Mask_B5';
% [~, Value.roiNames] = ReadTifFileNames(Value.roiMask);
% 
% Value.bgMask = 'I:\20190306_MoS2_1123_CH18SH\Mask_B5_bg';
% [~, Value.bgNames] = ReadTifFileNames(Value.bgMask);
% 
% Value.tifSeq = cell(size(Value.roiNames));
% tic
% for n = 1:length(Value.roiNames)
%     roiMask = ~imread(fullfile(Value.roiMask, Value.roiNames{n}));
%     bgMask = ~imread(fullfile(Value.bgMask, Value.bgNames{n}));
%     
%     check = xor(roiMask, bgMask);
%     if sum(check(:)) == 0
%         return
%     end
%     clear check
%     
%     for ii = Value.beginFrame:Value.endFrame
%         tif  = double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name))) - tif0;
%         Value.tifSeq{n, 1}((ii-Value.beginFrame+1), 1) = ROImean(tif, roiMask);
%         Value.tifSeq{n, 1}((ii-Value.beginFrame+1), 2) = ROImean(tif, bgMask);
%     end
%     
%     [f1, P1] = fft_P1(Value.tifSeq{n, 1}(:, 1), Fs);
%     figure('color', 'w');
%     plot(f1, P1);
%     hold on
%     [f2, P2] = fft_P1(Value.tifSeq{n, 1}(:, 2), Fs);
%     plot(f2, P2);
%     %     xlim([0, 20]); ylim([0, 50]);% set Fs = 5 Hz
%     xlim([0, 30]); ylim([0, 60]);
%     xlabel('f (Hz)','fontsize',10)
%     ylabel('|P(f)|','fontsize',10)
%     set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2)
%     set(gca, 'linewidth', 1.5)
%     title(n)
%     legend('ROI', 'Background')
%     hold off
%     
% end
% toc
% 
% %% Total sampleArea
% [cstrFilenames, cstrPathname] = uigetfile(...
%     {'*.*',  'All Files (*.*)';...
%     '*.zip',  'Zip-files (*.zip)';...
%     '*.roi',  'ROI (*.roi)'...
%     },'Pick a .roi imageJ file');
% [sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
% for m = 1:length(sROI)
%     [row, col] = ImageJroiLocation(sROI{m});
%     
%     Value.sampleArea = cell(size(Value.frame, 1), 2);
%     for ii = Value.beginFrame:Value.endFrame
%         tif  = double(imread(fullfile(Value.tifFile, Value.tifDir(ii).name))) - tif0;
% %         D3.sampleArea{(ii-D3.beginFrame) + 1, 1} = ii - D3.beginFrame + 1;
%         Value.sampleArea{(ii-Value.beginFrame) + 1, m} = tif((row(1):row(2)), (col(1):col(2)));
%     end
%     
% end
% 
%     
% %% Try adpmedian
% n = 20;
% 
% img = Value.sampleAreaFFT_f{n};
% % subplot(1, 2, 1);
% subplot(1, 3, 1);
% imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
% set(gca, 'CLim', [-60 60]);
% colormap jet
% % map = colormap('jet');
% % colorbar;
% % imshow(img, 'InitialMagnification', 'fit', 'colormap', hsv);
% impixelinfo
% 
% img1 = Value.sampleAreaFFT_f{n+1};
% % subplot(1, 2, 2);
% subplot(1, 3, 2);
% imshow(img1, 'DisplayRange', [], 'InitialMagnification', 'fit');
% set(gca, 'CLim', [-60 60]);
% colormap jet
% impixelinfo
% 
% subplot(1, 3, 3);
% imshow((img1 - img), 'DisplayRange', [], 'InitialMagnification', 'fit');
% set(gca, 'CLim', [-15 15]);
% colormap jet
% impixelinfo
% 
% %% Total sampleMask
% ii = 45;
% img = Value.sampleArea{ii, 1};
% figure('color', 'w');
% imshow(img)
% h_img = imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
% img1 = fft2(img);
% img2 = fftshift(img1);
% figure('color', 'w');
% imshow(img2, 'DisplayRange', [], 'InitialMagnification', 'fit');
% h = drawrectangle;
% Value.sampleMask = createMask(h, h_img);
% 
% 
% %% Total FFT
% 
% tic
% sampleArea = Value.sampleArea;
% sampleMask = Value.sampleMask;
% sampleAreaFFT = cell(size(sampleArea, 1), 1);
% parfor n = 1:size(sampleArea, 1)
%     sampleAreaFFT{n, 1} = FFTconvert(sampleArea{n, 1}, sampleMask);
% end
% Value.sampleAreaFFT =  sampleAreaFFT;
% clear  sampleAreaFFT  sampleArea
% toc
% 
% %% Figure
% 
% n = 91;
% img = Value.sampleArea{n, 2};
% subplot(2, 1, 1);
% imshow(img, 'DisplayRange', [], 'InitialMagnification', 'fit');
% set(gca, 'CLim', [-60 60]);
% colormap jet
% impixelinfo
% 
% m = n - 6;
% img1 = Value.sampleAreaFFT_f{m, 1};
% % img1 = B2.sampleAreaF1_3Df(:,:,n);
% % img1 = B2.sampleAreaF1_3Df2(:,:,n);
% subplot(2, 1, 2);
% imshow(img1,'DisplayRange', [],'InitialMagnification', 'fit')
% % imshow(img1, 'DisplayRange', [], 'InitialMagnification', 'fit');
% set(gca, 'CLim', [-40 40]);
% colormap jet
% % colorbar
% impixelinfo
% 
% %% Total adpmedian
% 
% Value.sampleAreaF1 = cell(size(Value.sampleArea));
% for n = 1:size(Value.sampleArea, 1)
%     img = Value.sampleArea{n, 2};
%     Value.sampleAreaF1{n, 1} = n;
%     Value.sampleAreaF1{n, 2} = adpmedian(img, 3);
% end
% 
% %% Try lowpass Analysis
% for n = 1:length(Value.tifSeq)
%     temp = Value.tifSeq{n, 1};
%     roi = lowp(temp(:, 1), 1, 36, 0.1, 20, Fs);
%     bg = lowp(temp(:, 2), 1, 36, 0.1, 20, Fs);
%     figure('color', 'w');
%     plot(Value.frame, roi, Value.frame, bg)
%     legend('roi', 'bg')
%     title(n)
%     xlabel('Frames')
%     ylabel('Intensity')
% end
% 
% %%
% U = lowp(I, 1, 4, 33*0.01, 20, 100);
% C = -intensity2current(U, 1601);
% plot(P, C);
% 
% %% Try Time-Frequency Analysis
% 
% temp = Value.tifSeq{n, 1};
% % roi = lowp(temp(:, 1), 4, 12, 0.1, 20, Fs);
% roi = lowp(temp(:, 1), 4, 16, 33*0.01, 20, Fs);
% roi = highpass(roi, 10, Fs);
% bg = lowp(temp(:, 2), 4, 16, 33*0.01, 20, Fs);
% % bg = lowp(temp(:, 2), 3, 7, 0.1, 20, Fs);
% bg = highpass(bg, 10, Fs);
% % figure('color', 'w');
% plot(x, roi, x, bg)
% % plot(x, temp(:, 1), x, temp(:, 2))
% legend('roi', 'bg')
% title(n)
% xlim([200 300])
% xlabel('Frames')
% ylabel('Intensity')
% figure('color', 'w');
% subplot(2,1,1)
% % spectrogram(temp(:, 1),[],[],[],Fs); title('roi');
% spectrogram(roi,[],[],[],Fs); title('roi');
% set(gca, 'CLim', [-40 40]);
% % colormap jet
% % figure('color', 'w');
% subplot(2,1,2)
% % spectrogram(temp(:, 2),[],[],[],Fs); title('roi');
% spectrogram(bg,[],[],[],Fs); title('bg');
% set(gca, 'CLim', [-40 40]);
% % colormap jet
% 
% %% Try one-dot Time-Frequency Analysis
% % from 2D to 3D
% 
% a = Value.sampleAreaF1{1,2};
% a(:,:,2) = Value.sampleAreaF1{2,2};
% 
% A = {zeros(10),ones(10)};
% reshape(cell2mat(A),[size(A{1}) numel(A)]);
% 
% %% Total Time-Frequency Analysis - fault
% 
% temp = Value.sampleAreaF1(:, 2);
% Value.sampleAreaF1_3D = reshape(cell2mat(temp),[size(temp{1}) numel(temp)]);
% clear temp
% 
% tic
% Value.sampleAreaF1_3Df = lowp_s(Value.sampleAreaF1_3D, Fs);
% toc
% 
% temp = Value.sampleAreaFFT;
% Value.sampleAreaF1_3D = reshape(cell2mat(temp),[size(temp{1}) numel(temp)]);
% clear temp
% Value.sampleAreaFFT1 = lowp_s(Value.sampleAreaF1_3D, Fs);
% 
% %% Total Time-Frequency Analysis
% 
% % temp = B2.sampleAreaFFT(6:1029);
% temp = Value.sampleAreaFFT;
% % sampleAreaFFT_3D = zeros([size(B2.sampleMask) 1024]);
% for n = 1:length(temp)
%     Value.sampleAreaFFT_3D(:, :, n) = temp{n, 1};
% end
% clear temp
% 
% row = size(Value.sampleAreaFFT_3D, 3);
% Value.sampleAreaFFT_3D_slide = zeros(row, 2);
% for n = 1:row
%     Value.sampleAreaFFT_3D_slide(n, 1) = ROImean(Value.sampleAreaFFT_3D(:, :, n), BW1);
%     Value.sampleAreaFFT_3D_slide(n, 2) = ROImean(Value.sampleAreaFFT_3D(:, :, n), BW2);
% end
% 
% tic
% Value.sampleAreaFFT_3Df = lowp_s(Value.sampleAreaFFT_3D, Fs);
% toc
% 
% tic
% temp1 = Value.sampleAreaFFT_3Df;
% sampleAreaFFT_f = cell(size(Value.sampleAreaFFT));
% parfor n = 1:row
%     sampleAreaFFT_f{n, 1} = temp1(:, :, n);
% end
% Value.sampleAreaFFT_f = sampleAreaFFT_f;
% clear temp1 sampleAreaFFT_f
% toc
% 
% %% Total differential intensity analysis 1
% BW1 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL1.tif');
% BW2 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL2.tif');
% BW3 = ~imread('I:\20190306_MoS2_1123_CH18SH\tifMin0_MaskL3.tif');
% row = size(Value.sampleAreaFFT_f, 1);
% Value.sampleAreaFFT_fL1L2 = zeros(row, 3);
% for n = 1:row
%     Value.sampleAreaFFT_fL1L2(n, 1) = ROImean(Value.sampleAreaFFT_f{n}, BW1);
%     Value.sampleAreaFFT_fL1L2(n, 2) = ROImean(Value.sampleAreaFFT_f{n}, BW2);
%     Value.sampleAreaFFT_fL1L2(n, 3) = ROImean(Value.sampleAreaFFT_f{n}, BW3);
% end
% L1_L2 = Value.sampleAreaFFT_fL1L2(:, 1) - Value.sampleAreaFFT_fL1L2(:, 2);
% L2_L3 = Value.sampleAreaFFT_fL1L2(:, 2) - Value.sampleAreaFFT_fL1L2(:, 3);
% L1_L3 = Value.sampleAreaFFT_fL1L2(:, 1) - Value.sampleAreaFFT_fL1L2(:, 3);
% x1 = ceil(Value.begin.pike+(Value.beginFrame/100*10000)-1);
% x2 = ceil(Value.begin.pike+(Value.beginFrame/100*10000)-1+1023/100*10000);
% Value.t1 = (x1/10000:0.0001:x2/10000)';
% Value.t1_y = Value.data(x1:x2, 2);
% x3 = (Value.beginFrame:1:Value.endFrame)';
% Value.t2 = zeros(1024, 1);
% for ii = 1:1024
%     Value.t2(ii) = (Value.begin.pike+x3(ii)/106*10000-1)/10000;
% end
% figure('color', 'w');
% subplot(2, 1, 1);
% plot(Value.t1, 10*Value.t1_y, '-', Value.t2, L1_L2, Value.t2, L2_L3, Value.t2, L1_L3)
% subplot(2, 1, 2);
% plot(Value.t1, 10*Value.t1_y, '-', Value.t2, Value.sampleAreaFFT_fL1L2(:, 1), Value.t2, Value.sampleAreaFFT_fL1L2(:, 2), Value.t2, Value.sampleAreaFFT_fL1L2(:, 3))
% 
% %% Try differential intensity analysis 2
% Value.sampleAreaFFT_f_diff = cell(size(Value.sampleAreaFFT_f, 1)-1,1);
% for ii = 1:(1024-1)
%     Value.sampleAreaFFT_f_diff{ii} = Value.sampleAreaFFT_f{ii+1} - Value.sampleAreaFFT_f{ii};
% end
% 
% %%
% Fs = 100;
% N = 1024;
% n = 0:(N-1) ;
% y = fft(x, N);
% 
% %% Total get Phase difference of Voltage vs time and Intensity vs time
% n = 1024;
% skip = length(Value.t1)/n;
% Value.t1_y0 = Value.t1_y(1:skip:end);
% Value.PhDiff = phdiffmeasure(Value.t1_y0(1:1024), L1_L2);
% 
% %% Try phase difference 1
% [Amp2, Ph2] = freqrespmeasure(B2.t1_y0(1:1024), B2.L1_L2);
% [Value.Amp, Value.Ph] = freqrespmeasure(t1_y0(1:1024), Value.L1_L2);
% 
% 
% %% Total 3D phase difference
% tic
% temp = Value.sampleArea;
% % sampleAreaFFT_3D = zeros([size(B2.sampleMask) 1024]);
% for n = 1:length(temp)
%     Value.sampleArea_3D(:, :, n) = temp{n, 1};
% end
% clear temp
% 
% output = freqrespmeasure_s(Value.sampleArea_3D, Value.t1_y0);
% % Ma = max(output(:));
% % Mi = min(output(:));
% % OneMi = ones(size(B5.sampleMask))*Mi;
% output = output + (output<0)*2*pi - (output>2*pi)*2*pi;
% figure('color', 'w');
% imshow(output, 'DisplayRange', [], 'InitialMagnification', 'fit');
% % set(gca, 'CLim', [Mi Ma]);
% colormap jet
% colorbar
% impixelinfo
% Value.phaseDiff = output;
% toc
% %% Try Ampliture 2
% tic
% 
% % temp = D3.sampleAreaFFT;
% % for n = 1:length(temp)
% %     D3.sampleAreaFFT_3D(:, :, n) = temp{n, 1};
% % end
% % clear temp
% 
% % setFrequency = 5;
% output = extremum_Amp(Value.sampleArea_3D);
% figure('color', 'w');
% imshow(output, 'DisplayRange', [], 'InitialMagnification', 'fit');
% colormap jet
% colorbar
% impixelinfo
% Value.Ampliture = output;
% toc
% %%
% figure('color','w');
% subplot(5, 1, 1); hist(DphaseDiff.D3_5)
% subplot(5, 1, 2); hist(DphaseDiff.D4_10)
% subplot(5, 1, 3); hist(DphaseDiff.D5_15)
% subplot(5, 1, 4); hist(DphaseDiff.D1_20)
% subplot(5, 1, 5); hist(DphaseDiff.D2_25)
% 
% %%
% toc
% disp(['Running time: ',num2str(toc)]);
% 
% MailToMe('nona1588@outlook.com');
% 
% 
% %% Total the experiment of Ru(III)
% tifFile = 'E:\20181116_MoS2_CH18-Au\D3_Z1_Ru_PBS_0 -0-4V_0-1VpS_2c_HMMT_200fps';
% [Value.tifFolder, Value.tifNames] = ReadTifFileNames(tifFile);
% Value.tif0 = double(imread(fullfile(Value.tifFolder, Value.tifNames{1})));
% [cstrFilenames, cstrPathname] = uigetfile(...
%     {'*.*',  'All Files (*.*)';...
%     '*.zip',  'Zip-files (*.zip)';...
%     '*.roi',  'ROI (*.roi)'...
%     },'Pick a .roi imageJ file');
% [sROI] = ReadImageJROI(fullfile(cstrPathname, cstrFilenames));
% m = 2;
% Value.beginFrame = 417;
% Value.endFrame = 2017;
% % C1.sampleArea = cell(C1.endFrame - C1.beginFrame + 1,1);
% Value.sampleArea = cell(1600/20+1,2);
% [row, col] = ImageJroiLocation(sROI{m});
% n = 1;
% for ii = Value.beginFrame:20:Value.endFrame
%     tif  = double(imread(fullfile(Value.tifFolder, Value.tifNames{ii}))) - Value.tif0;
%     Value.sampleArea{n, 1} = n;
%     Value.sampleArea{n, 2} = tif((row(1):row(2)), (col(1):col(2)));
%     n = n + 1;
% end
% %%
% D3_C1.sampleArea = cell(81, 1);
% for ii = 1:81
%     D3_C1.sampleArea{ii} = Value.sampleArea{ii, 2} - C1.sampleArea{ii, 2};
% end