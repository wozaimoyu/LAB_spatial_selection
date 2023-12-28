% Построение рис. 1, 2 и 4

clear all; close all; clc;
mkdir pics

% Цветовая гамма рисунков
% 1 - цветные рисунки, 0 - ч\б рисунки
colorPlt = 0;
% автоматическое сохранение рисунков
saveFig = 1;



stdCoordsArr = [0, 5];  % СКО оценки координат UE по осям x, y, z в метрах
Nel = 8;
% выбор типа антенной решетки
% 1 - планарная АР
% 2 - линейная АР
% 3 - круговая АР (поддержаны только режимы ДН antPattCntrl=0,1,2)
antType = 1;
antPattCntrlArr = [0, 1, 2];

pltStyle = ["k-"; "k--"];
ZZa = [];
for s=1:length(stdCoordsArr)
    ZZ = [];
    stdCoords = stdCoordsArr(s);
    for aa=antPattCntrlArr
        antPattCntrl = aa;
        if aa ~= 2
        rng('default')
        else
        rng(s)
        end
        lab_spatialfilt_base;
        ZZ = [ZZ; Z];
    end
    ZZa = [ZZa; Z];
    figure;
    hold on
    for ii=1:2
        if colorPlt == 1
            plot(X, ZZ(ii,:))
        else
            plot(X, ZZ(ii,:), pltStyle(ii))
        end
    end
    ax = gca;
    if colorPlt == 0
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    ax.GridColor = [0 0 0];
    end
    grid on;     xlabel('x, м'); ylabel('SIR_{avg}, дБ');
    legend('max ДН', 'max+null ДН')
    % автоматическое сохранение рисунков
    if saveFig  == 1
    exportgraphics(ax, sprintf('pics/рис%i.tif', s),'Resolution',300);
    end

end
figure;
hold on
for ii=1:2
    if colorPlt == 1
        plot(X, ZZa(ii,:))
    else
        plot(X, ZZa(ii,:), pltStyle(ii))
    end
end
ax = gca;
if colorPlt == 0
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.GridColor = [0 0 0];
end
grid on; xlabel('x, м'); ylabel('SIR_{avg}, дБ');
legend(sprintf("СКО=%iм",stdCoordsArr(1)),...
       sprintf("СКО=%iм",stdCoordsArr(2)));
% автоматическое сохранение рисунков
if saveFig == 1
exportgraphics(ax, 'pics/рис4.tif','Resolution',300);
end

