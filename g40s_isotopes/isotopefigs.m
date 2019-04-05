%% Set seawater conditions
psal = 35;
d13dic = 2; % permille
talk = 2325; % umol / kg
si = 0; % umol / kg
phos = 0; % umol / kg

%% Draw figure
figure(1); clf
printsetup(gcf, [10 10])
subplot(2, 1, 1); hold on
    temps = (0:35)';
    rs = [5 10 20] * 1e-6;
    rclr = {'k' 'r' 'b' 'm'};
    pco2s = [400 500];
    psty = {'-' '--'};
    d13pocs = NaN(numel(temps), numel(pco2s), numel(rs));
    eps = NaN(numel(temps), numel(pco2s), numel(rs));
    for R = 1:numel(rs)
        for P = 1:numel(pco2s)
            [d13pocs(:, P, R), eps(:, P, R)] = ...
                rau1996(rs(R), temps, psal, talk, pco2s(P), ...
                d13dic, si, phos);
            plot(temps, eps(:, P, R), ...
                'color', rclr{R}, ...
                'linestyle', psty{P})       
        end % for P
    end % for R
    xlabel(['Temperature / ' degrees 'C'])
    ylabel([grk('e') '_p / ' permille])
    setaxes(gca, 8)
    xlim([0 35])
    ylim([6 24])
    set(gca, ...
        'box', 'on', ...
        'xtick', 0:7:35, ...
        'ytick', 6:6:24)
    text(0, 1.1, '(a)', ...
        'FontName', 'arial', ...
        'FontSize', 8, ...
        'Units', 'normalized')
    ax1 = gca;
subplot(2, 1, 2); hold on
    d13pocs2 = NaN(numel(temps), numel(rs));
    co2s1 = CO2SYSv2(talk, pco2s(1), 1, 4, psal, temps, temps, 0, 0, ...
        si, phos, 3, 10, 3);
    co2s2 = CO2SYSv2(talk, pco2s(2), 1, 4, psal, temps, temps, 0, 0, ...
        si, phos, 3, 10, 3);
    dic1 = co2s1(:, 2);
    dic2 = co2s2(:, 2);
    deldic = dic2 - dic1;
    delRC = -0.016; % permille / (umol / kg)
    deld13dic = deldic * delRC;
    fleg = cell(numel(pco2s) * numel(rs), 1);
    LC = 1;
    for R = 1:numel(rs)
        % First draw fractionation-only change
        for P = 1:numel(pco2s)
            plot(temps, d13pocs(:, P, R), ...
                'color', rclr{R}, ...
                'linestyle', psty{P})
            fleg{LC} = [num2str(pco2s(P)) ' ppm, ' num2str(rs(R) * 1e6) ...
                ' ' grk('m') 'm'];
            LC = LC + 1;
        end % for P
        % Now add changes in everything else
        d13pocs2(:, R) = ...
            rau1996(rs(R), temps, psal, talk, pco2s(2), ...
            d13dic + deld13dic, si, phos);
        plot(temps, d13pocs2(:, R), ...
            'color', rclr{R}, ...
            'linestyle', ':')
        fleg{LC} = 'As above + S.E.';
        LC = LC + 1;
    end % for R
    xlabel(['Temperature / ' degrees 'C'])
    ylabel([grk('d') '^{13}C_{POC} / ' permille])
    fl = legend(fleg, 'location', 'eastoutside');
    fl.FontName = 'arial';
    fl.FontSize = 8;
    setaxes(gca, 8)
    xlim([0 35])
    ylim([-33 -13])
    set(gca, ...
        'box', 'on', ...
        'XTick', 0:7:35, ...
        'YTick', -33:5:-10)
    text(0, 1.1, '(b)', ...
        'FontName', 'arial', ...
        'FontSize', 8, ...
        'Units', 'normalized')
    ax2 = gca;
ax1.Position = [0.12 0.6 0.5 0.34];
ax2.Position = [0.12 0.1 0.5 0.34];
fl.Position = [0.7 0.4 0.25 0.2];
% print('-r300', ...
%     'E:\Dropbox\Write-ups\GEOTRACES 40S\Tuerena - isotopes\rau1996', ...
%     '-dpng')

%% Get Suess effect numbers
dfrac = squeeze(diff(d13pocs, 1, 2));
dboth = d13pocs2 - squeeze(d13pocs(:, 1, :));
dsuess = dboth - dfrac;
pctsuess = 100 * dsuess ./ dboth;

%% Climate sensitivity calculation
cspco2 = 100:1000;
cstemp = logn(cspco2, 2) * 1.5;
figure(2); clf
    plot(cspco2, cstemp)
    grid on
dt54 = (logn(500, 2) - logn(400, 2)) * [1.5 4.5];
[d5d13pocs, d5eps] = rau1996(10e-6, 17, psal, talk, 500, ...
    d13dic, si, phos);
[d4d13pocs, d4eps] = rau1996(10e-6, 17, psal, talk, 400, ...
    d13dic, si, phos);
[d4td13pocs, d4teps] = rau1996(10e-6, 17+dt54', psal, talk, 400, ...
    d13dic, si, phos);
d54d13_pco2 = d5d13pocs - d4d13pocs;
d54d13_temp = d4td13pocs - d4d13pocs;
