%% -- raw figures
Potential = sheet(:, 10);
Current = zeros(size(Yfit(2:end, :)));
for n = 1:col
    Current(:, n) = intensity2current(Yfit(:, n), 1601); % 1601 is the number of Yfit(:, n)
end

figure('color','w');
for n = 1:col
    plot(Potential(4:end), -Current(3:end, n)); % get fitted lines
    hold on
end

xlabel('Potential/V','fontsize',10);
ylabel('I-Current','fontsize',10);
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
set(gca, 'linewidth', 1.5)

%% -- signal frequency spectrum
% Fs = 100; % REAL hamamatsu camera in CV
% Fs = 800; % hamamatsu camera in CV
Fs = 100; % normal sampling of Pike camera
% Fs = 200; % hamamatsu camera in i-t
fX1 = intensity;
% fX2 = Y(:, 2);
% fX3 = Y(:, 3);
% Y2 = fft(X2);
% L = length(X2);
% P2 = abs(Y2/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
[f1, P1] = fft_P1(fX1, Fs);
% [f2, P2] = fft_P1(fX2, Fs);
% [f3, P3] = fft_P1(fX3, Fs);
figure('color','w');
% subplot(3,1,1);
plot(f1, P1);
xlim([1, 50]);
% ylim([0, 6]);
% title('Single-Sided Amplitude Spectrum of Graphene blinking (HMMT)')
% legend('Background: Au');
% hold on
% subplot(3,1,2);
% plot(f2, P2);
% xlim([0, 53]); ylim([0, 3]);
% legend('Graphene, near the edge');
% subplot(3,1,3);
% plot(f3, P3);
% xlim([0, 53]); ylim([0, 3]);
% legend('on Graphene');
% title('Single-Sided Amplitude Spectrum of Graphene blinking')
xlabel('f (Hz)','fontsize',10)
ylabel('|P(f)|','fontsize',10)
% set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
% set(gca, 'linewidth', 1.5)

%%
color = [...
    'r',... % Red 1
    'g',... % Green 2
    'b',... % Blue 3
    'c',... % Cyan 4
    'm',... % Magenta 5
    'y',... % Yellow 6
    'k',... % Black 7
    'w',... % White 8
    ];
roiLabel = ['1', '2', '3', '4',...
    '5', '6', '7'];

%% -- filtered 1-D digital i-t signal
filtY = zeros(size(Y));
for n = 1:col
    filtY(:, n) = lowp(Y(:, n), 1, 4, 0.1, 20, Fs);
end

figure('color','w');
for n = 1:col
    plot(X, filtY(:, n), color(n));
    hold on
end

legend
xlabel('Frame');
ylabel('¦¤Intensity');

%% -- filtered 1-D digital CV signal
filtCurrent = zeros(size(Current(3:end, :)));
for n = 1:col
    filtCurrent(:, n) = lowp(Current(3:end, n), 4, 30, 0.1, 20, Fs);
end

figure('color','w');
for n = 1:col
    % scatter(Potential(4:end), -filtCurrent(:, n));
    %     figure('color','w');
    plot(Potential(4:end), -filtCurrent(:, n), color(n));
    %     legend(roiLabel(n)); % plot respectively
    %     xlabel('Potential/V');
    %     ylabel('¦¤Current');
    hold on
end

legend
xlabel('Potential/V','fontsize',10);
ylabel('I-Current','fontsize',10);
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 2);
set(gca, 'linewidth', 1.5)
%%
plot(X, intensity); 
hold on
intensity2 = lowp(intensity, 1, 36, 0.1, 20, 100);
plot(X, intensity2); 
ylim([1.905*10^8, 1.935*10^8]);
hold off

%% for TaS2's CV
filtY2= lowp(Y(:,2), 1, 33, 0.1, 18, 106);
filtY2=smooth(filtY2, 5);
subplot(1,2,1);
% plot(X(849:1697), filtY2(849:1697));
plot(X(1:849), filtY2(1:849));
xlim([-0.8,0]);
xlabel('Potential/V','fontsize',10);
ylabel('¦¤intensity/IU','fontsize',10);
% hold on
% diffY2 = diff(filtY2);
diffY2 = lowp(diff(filtY2), 2, 10, 0.1, 20, 106);
subplot(1,2,2);
% plot(X(849:1696), -10*diffY2(849:1696));
plot(X(1:848), -10*diffY2(1:848));
xlim([-0.8,0]);
xlabel('Potential/V','fontsize',10);
ylabel('First-order derivative of SPR intensity/(IU/s)','fontsize',10);
% filtC2= lowp(C2(3:end), 4, 23, 0.1, 20, 800);
% plot(Potential(4:801), -filtC2(1:798));
% hold off
%%
for n = 1:col
    figure('color', 'w');
    plot(X(2:end), diffY(:, n));
    hold on
    plot(X(2:end), f_dY(:, n));
    hold off
end
