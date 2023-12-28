% Построение рис. 3, 14, 15

clear all; close all; clc;
mkdir pics

% Цветовая гамма рисунков
% 1 - цветные рисунки, 0 - ч\б рисунки
colorPlt = 0;
% автоматическое сохранение рисунков
saveFig = 1;



stdCoordsArr = [2, 5, 10, 15];  % СКО оценки координат UE по осям x, y, z в метрах
NellArr = 3:3:30;
% выбор типа антенной решетки
% 1 - планарная АР
% 2 - линейная АР
% 3 - круговая АР (поддержаны только режимы ДН antPattCntrl=0,1,2)
antType = 1;


pltStyle_1 = ["ko-"; "k+-"; "k*-"; "k^-"];
pltStyle_2 = ["ko--"; "k+--"; "k*--"; "k^--"];
picName = ["рис3", "рис15", "рис14"];

% массив с номерами режимов ДО, адаптивная схема не строится, тк не зависит
% от СКО
antPattCntrlArr = [0, 1, 3, 3];
win_typeArr = [0, 0, 0, 2];
antPattCntrlCmt = ["Управление максимумом ДН", ...
                   "Управление максимумом и нулем ДН", ...
                   "Управление шириной ДН Гаусса", ...
                   "Управление шириной прямоугольной ДН"];
antTypeCmt = ["Планарная АР", "Линейная АР", "Круговая АР"];
Zaa = [];
for aa=1:length(antPattCntrlArr)
Zst = [];
antPattCntrl = antPattCntrlArr(aa);
win_type = win_typeArr(aa);
for s=stdCoordsArr
    Zn = [];
    stdCoords = s;
    for nel=NellArr
        Zsc = [];
        Nel = nel;
        switch win_type
            case 0
                sc = 3.79*exp(-0.126*stdCoords);
            case 1
                sc = 0.1;
            case 2
                sc = -0.2*stdCoords + 4.65;
        end
        rng('default')
        fprintf('Alg %i std=%.2f, Nel=%i, Scl=%.2f\n', antPattCntrl, s, nel, sc);
        lab_spatialfilt_base;
        Zsc = [Zsc, mean(Z(:))];
        Zn = [Zn;Zsc];
    end
    Zst{end+1} = Zn;
end
Zaa{end+1} = Zst;
end

%%
ZZmS_1 = Zaa{1};
for aa=2:length(antPattCntrlArr)
    ZZmS = Zaa{aa};
    figure; hold on; grid on; axis tight;
    for kk=1:size(ZZmS,2)
        if colorPlt == 1
            plot(NellArr, ZZmS_1{kk}, '-');
        else
            plot(NellArr, ZZmS_1{kk}, pltStyle_1(kk));
        end
    end
    set(gca,'ColorOrderIndex',1)
    for kk=1:size(ZZmS,2)
        if colorPlt == 1
            plot(NellArr, ZZmS{kk}, '--');
        else
            plot(NellArr, ZZmS{kk}, pltStyle_2(kk));
        end
    end
    ax = gca;
    if colorPlt == 0
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    ax.GridColor = [0 0 0];
    end
    xlabel('N');
    ylabel('SIR_{avg}, дБ');
    legFrst = ['СКО=',num2str(stdCoordsArr(1)), 'м'];
    if antPattCntrlArr(aa) == 3
    legend([legFrst, strcat(string(stdCoordsArr(2:end)),' м')],'NumColumns',3, 'Location', 'southeast');
    else
    legend([legFrst, strcat(string(stdCoordsArr(2:end)),' м')],'NumColumns',3);
    end
    title([antPattCntrlCmt(aa), antTypeCmt(antType)]);
    if saveFig  == 1
    exportgraphics(ax, strcat('pics/', picName(aa-1), 'а.tif'),'Resolution',300);
    end
    
    figure; hold on; grid on; axis tight;
    for kk=1:size(ZZmS,2)
        if colorPlt == 1
            plot(NellArr, ZZmS{kk}-ZZmS_1{kk}, '-');
        else
            plot(NellArr, ZZmS{kk}-ZZmS_1{kk}, pltStyle_1(kk));
        end
    end
    ax = gca;
    if colorPlt == 0
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    ax.GridColor = [0 0 0];
    end
    xlabel('N');
    ylabel('\DeltaSIR_{avg}, дБ');
    legFrst = ['СКО=',num2str(stdCoordsArr(1)), 'м'];
    if antPattCntrlArr(aa) == 3
    legend([legFrst, strcat(string(stdCoordsArr(2:end)),' м')],'NumColumns',3, 'Location', 'southeast');
    else
    legend([legFrst, strcat(string(stdCoordsArr(2:end)),' м')],'NumColumns',3);
    end
    title([antPattCntrlCmt(aa), antTypeCmt(antType)]);
    if saveFig  == 1
    exportgraphics(ax, strcat('pics/', picName(aa-1), 'б.tif'),'Resolution',300);
    end
end
