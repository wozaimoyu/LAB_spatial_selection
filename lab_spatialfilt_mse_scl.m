% Построение рис. 12, 13

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


picName = ["рис12", "рис13"];


% выбор формы ДН для алгоритма управления шириной ДН (antPattCntrl = 3)
% 0 - окно Гаусса
% 1 - окно приподнятого косинуса
% 2 - прямоугольное окно
win_typeArr = [2, 0];
% массив с номерами режимов ДО
antPattCntrlArr = [3, 3];
antPattCntrlCmt = ["Управление максимумом ДН", ...
                   "Управление максимумом и нулем ДН", ...
                   "Адаптивное управление ДН", ...
                   "Управление шириной ДН Гаусса", ...
                   "Управление шириной прямоугольной ДН"];
antTypeCmt = ["Планарная АР", "Линейная АР", "Круговая АР"];
Zaa = [];
for aa=1:length(antPattCntrlArr)
Zst = [];
antPattCntrl = antPattCntrlArr(aa);
win_type = win_typeArr(aa);
switch win_type
    case 0
        sclArr = 0.1:0.5:5.5; % окно Гаусса
    case 1
        sclArr = 0:0.1:0.9; % окно raised cosine (не исследовалось)
    case 2
        sclArr = 0.5:0.5:5.5; % прямоугольное окно
end
for s=stdCoordsArr
    Zn = [];
    stdCoords = s;
    for nel=NellArr
        Zsc = [];
        Nel = nel;
        for sc=sclArr
            rng('default')
            fprintf('Win %i std=%.2f, Nel=%i, Scl=%.2f\n', win_type, s, nel, sc);
            lab_spatialfilt_base;
            Zsc = [Zsc, mean(Z(:))];
        end
        Zn = [Zn;Zsc];
    end
    Zst{end+1} = Zn;
end
Zaa{end+1} = Zst;
end

%%
for aa=1:length(antPattCntrlArr)
figure;
ZZmS = Zaa{aa};
[X,Y] = meshgrid(NellArr, sclArr);
for kk=1:size(ZZmS,2)
    subplot(2,2,kk)
    if colorPlt == 1
        surf(X, Y, ZZmS{kk}.', 'FaceColor', 'interp', 'EdgeColor','none');
        c1 = colorbar; c1.Label.String = 'SIR_{avg}, дБ';
    else
        mm = ceil(max(ZZmS{kk}(:)));
        mn = ceil(min(ZZmS{kk}(:)));
        lev = ceil(mn:(mm-mn)/8:mm);
        [~,hh] = contour(X, Y, (ZZmS{kk}.'),...
            lev, 'k-', 'ShowText', 'on');
        hh.LabelSpacing = 250;
        ax = gca;
        ax.XColor = [0 0 0];
        ax.YColor = [0 0 0];
        ax.GridColor = [0 0 0];
    end
    grid on; xlabel('N'); ylabel('s'); view([0, 90]);    
    axis tight
    title(sprintf('СКО=%d м',stdCoordsArr(kk)));
end
if win_typeArr(aa) == 0
    sgtitle(sprintf('ДН Гаусса %s', antTypeCmt(antType)));
elseif win_typeArr(aa) == 2
    sgtitle(sprintf('Прямоугольная ДН %s', antTypeCmt(antType)));
end
fg = gcf;
if saveFig  == 1
exportgraphics(fg, strcat('pics/', picName(aa), '.tif'),'Resolution',300);
end
end
