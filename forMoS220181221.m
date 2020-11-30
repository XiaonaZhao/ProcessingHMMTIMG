%% for MoS220181221
% add a New Variable, from
% 'E:\MoS2\20181221_MoS2_1219_BareAu\Result\A4_A3.opju', the Book2 -
% intensity

current = -intensity2current(intensity, 401)/100; % the Book4


%%
load('E:\MoS2\20181130_MoS2_CH18-Au-unclear-normal\DifferentConcentration.mat')

% figure('color', 'w');
% plot((1:length(A0))', A0)
% hold on
plot((1:length(B1))', lowp(B1, 1, 15, 0.001, 20, 100) )
hold on
plot((1:length(C5))', lowp(C5, 1, 15, 0.001, 20, 100))
plot((1:length(D10))', lowp(D10, 1, 15, 0.001, 20, 100))
% plot((1:length(E20))', E20)
hold off

%%
mM1 = 2*(B1(1004:1404, 1) - B1(1004)*ones(401, 1));
mM5 = 2*(C5(1027:1427, 1) - C5(1027)*ones(401, 1));
mM10 = 2*(D10(1012:1412, 1) - D10(1012)*ones(401, 1));
figure('color', 'w');
plot((1:length(mM1))', mM1)
hold on
plot((1:length(mM5))', mM5)
plot((1:length(mM10))', mM10)
xlim([0 400])
hold off


