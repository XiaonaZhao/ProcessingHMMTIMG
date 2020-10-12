% %% Figure 3 multilayer
% load('G:\MoS2\MoS2_0802\_Result_f\matlab_A3A4A5.mat')
% %%
% L = 850;
% 
% [row, col] = ImageJroiLocation(Value(1).sROI);
% tif = double(imread(fullfile(Value(1).tifPath, Value(1).tifName{1})));
% tif1 = tif(row(1):row(2), col(1):col(2));
% 
% exp = cell(L, 1);
% for ii = 1:L
%     tif = double(imread(fullfile(Value(1).tifPath, Value(1).tifName{ii+Value(1).begin-1})));
%     exp{ii, 1} = tif(row(1):row(2), col(1):col(2)) - tif1;
% end
load('G:\MoS2\MoS2_0802\_Result_f\A4_concentration.mat') % A4
% load('G:\MoS2\MoS2_0802\_Result_f\A3_concentration.mat') % A3 dizzy
%%
figure('color', 'w');

%% Cutscale Vision
savepath = 'G:\MoS2\MoS2_0802\_Result_f\TIF_A4_fire_RuII\';
for ii = 400:800
% ii = 600;
    tif_ii = exp{ii, 1} - exp{801, 1};
%     tif_ii = tif_ii(row(1):row(2), col(1):col(2))/2;
    
    if ii > 600
        Voltage = -0.3 + ((ii-600)/100)*0.1;
    else
        Voltage = -0.1 - ((ii-400)/100)*0.1;
    end
    
    redCon_ii = (-tif_ii)*10/(-(filterCurve(600)-filterCurve(400)));
    localSums = imboxfilt(redCon_ii, 13);
    
    localSums(localSums > 10) = 10;
    localSums(localSums < 0) = 0;
    
    imshow(localSums, 'DisplayRange',[], 'InitialMagnification', 'fit');
    title(['t = ' num2str((ii-400)/100), ' s, Voltage = ' num2str(Voltage), ' V'], 'FontSize', 14, 'FontWeight', 'bold');
    colormap fire % parula / turbo
    colorbar;

    set(gca, 'CLim', [0 10]);
    h = colorbar;
    set(h, 'FontSize', 14, 'FontWeight', 'bold');
    set(get(h,'title'),'string','[Ru(NH_3)_6]Cl_3^2^+ (mM)', 'FontSize', 12);

    pause(0.05);
    saveas(gcf,[savepath, 'A4_' num2str(ii, '%04d'), '.tif']);
    
end
