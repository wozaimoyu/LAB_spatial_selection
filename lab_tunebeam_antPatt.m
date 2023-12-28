% Построение рис. 5, 6, 7

clear all; close all; clc;
mkdir pics

% Цветовая гамма рисунков
% 1 - цветные рисунки, 0 - ч\б рисунки
colorPlt = 0;
% автоматическое сохранение рисунков
saveFig = 1;

picName = ["рис5", "рис6", "рис7"];


c = physconst('LightSpeed');
f = 30e9;       % несущая в диапазоне ММВ, Гц
lamb = c/f;     % длина волны, м
da = 0.5*c/f;   % расстояние между элементами АР

BW = 50; % ширина главно лепестка ДН, градусы
Nel = 20; % число АЭ в одном измерении
% выбор типа антенной решетки:
% 1 - планарная АР
% 2 - линейная АР
% 3 - круговая АР (не поддерживается)
antType = 1;
antElPos = createAnt(antType, Nel, da); % формирование АР

figNumber = 1;
for j=[0,1]
% расчет коэфф. для прямоугольной ДН
[wr, azAngP, antPattPr] = beamshapingWeight(2, BW, 0, Nel, 1, j);
if (antType == 1)
    % расчет вектора весовых коэфф. планарной АР для вертикальных АЭ
    wr = repmat(wr, Nel, 1)/Nel;
    wr = wr(:);
end
gr = zeros(1, length(azAngP));
for i=1:length(azAngP)
    gr(i) = getAntPatternG(antElPos, f, azAngP(i), 0, wr, 0);
end

figure;
if colorPlt == 1
    plot(azAngP, antPattPr);
    hold on
    plot(azAngP, gr);
else
    plot(azAngP, antPattPr, 'k-');
    hold on
    plot(azAngP, gr, 'k--');
    ax = gca;
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    ax.GridColor = [0 0 0];
end
grid on
xlabel('\phi, градус');
ylabel('|A(\phi)|');
legend('задано', 'синтез')
if saveFig  == 1
exportgraphics(figure(figNumber), strcat('pics/', picName(figNumber), '.tif'),'Resolution',300);
end
figNumber = figNumber + 1;
end

% расчет коэфф. для ДН Гаусса
[wg, azAngP, antPattPg] = beamshapingWeight(0, BW, 0, Nel, 1);
if (antType == 1)
    % расчет вектора весовых коэфф. планарной АР для вертикальных АЭ
    wg = repmat(wg, Nel, 1)/Nel;
    wg = wg(:);
end
gg = zeros(1, length(azAngP));
for i=1:length(azAngP)
    gg(i) = getAntPatternG(antElPos, f, azAngP(i), 0, wg, 0);
end

figure;
if colorPlt == 1
    plot(azAngP, antPattPg, 'LineWidth', 2);
    hold on
    plot(azAngP, gg, '--', 'LineWidth', 2);
else
    plot(azAngP, antPattPg, 'k-');
    hold on
    plot(azAngP, gg, 'k--', 'LineWidth', 3, 'MarkerIndices', 1:1e2:length(gg));
    ax = gca;
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    ax.GridColor = [0 0 0];
end
grid on
xlabel('\phi, градус');
ylabel('|A(\phi)|');
legend('задано', 'синтез')
if saveFig  == 1
exportgraphics(figure(figNumber), strcat('pics/', picName(figNumber), '.tif'),'Resolution',300);
end
