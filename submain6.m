
figure('color', 'w');
plot(potential, intensity)
xlabel('Potential (V vs. Ag/AgCl)')
ylabel('\Delta SPR intensity (a.u.)')
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 1);
set(gca, 'linewidth', 1, 'FontSize', 16)

%%
kB = 1.380649e-23;
T = 298.1;
electron = 1.6e-19;

x = (0:135)'/13.6;
kBT2elnR2O = zeros(136, 1);
for ii = 1:136
    kBT2elnR2O(ii) = kB*T/electron*log(x(ii)/(10-x(ii)));
end

X = potential(101:236);
figure('color', 'w');
%%
plot(-kBT2elnR2O, intensity(101:236))

xlabel('-RT/F·ln([Ru^2^+]/[Ru^3^+]) (V)')
ylabel('\Delta SPR intensity (a.u.)')
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 1);
set(gca, 'linewidth', 1, 'FontSize', 16)
%%
figure('color', 'w');
yyaxis left
plot(-kBT2elnR2O, intensity(101:236))
xlabel('-RT/F·ln([Ru^2^+]/[Ru^3^+]) (V)')
ylabel('\Delta SPR intensity (a.u.)')
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 1);
set(gca, 'linewidth', 1, 'FontSize', 16)
yyaxis right
plot(-kBT2elnR2O, x)
ylabel('[Ru^2^+] (mM)')
%%
RuIII = 10*ones(length(x), 1) - x;
figure('color', 'w');
yyaxis left
plot(-kBT2elnR2O, intensity(101:236))
xlabel('-RT/F·ln([Ru^2^+]/[Ru^3^+]) (V)')
ylabel('\Delta SPR intensity (a.u.)')
set(findobj(get(gca, 'Children'), 'LineWidth',0.5), 'LineWidth', 1);
set(gca, 'linewidth', 1, 'FontSize', 16)
yyaxis right
plot(-kBT2elnR2O, RuIII)
ylabel('[Ru^3^+] (mM)')