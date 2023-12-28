% Построение рис. 10, 11

clear all; close all; clc;
mkdir pics

% Цветовая гамма рисунков
% 1 - цветные рисунки, 0 - ч\б рисунки
colorPlt = 0;
% автоматическое сохранение рисунков
saveFig = 1;

picName = ["рис10", "рис11"];





stdCoordsArr = [2, 5, 10, 15];  % СКО оценки координат UE по осям x, y, z в метрах
NellArr = 3:3:30; % число АЭ в одном измерении
% выбор типа антенной решетки:
% 1 - планарная АР
% 2 - линейная АР
% 3 - круговая АР (не поддерживается)
antType = 1;





% выбор формы ДН для алгоритма управления шириной ДН (antPattCntrl = 3)
% 0 - окно Гаусса
% 1 - окно приподнятого косинуса
% 2 - прямоугольное окно
win_typeArr = [2, 0];
antTypeCmt = ["Планарная АР", "Линейная АР", "Круговая АР"];

for ww=1:length(win_typeArr)
win_type = win_typeArr(ww);
switch win_type
    case 0
        sclArr = 0.2:0.4:3; % окно Гаусса
    case 1
        sclArr = 0:0.1:0.9; % окно raised cosine
    case 2
        sclArr = 1:0.2:3; % прямоугольное окно
end

c = physconst('LightSpeed');
f = 30e9;       % несущая в диапазоне ММВ, Гц
lamb = c/f;     % длина волны, м
da = 0.5*c/f;   % расстояние между элементами АР

Z = [];

for s=stdCoordsArr
Zn = [];
for nel=NellArr
Zsc = [];
for sc=sclArr
rng('default')
    
stdCoords = s; % СКО оценки координат x, y, z

Nel = nel;                              % число АЭ в одном измерении
antElPos = createAnt(antType, Nel, da); % формирование АР
NelFull = size(antElPos, 1);

N = 100;  % число точек расчета
HgNB = 0; % высота подвеса АР gNB
Due = 50; % расстояние между gNB и UE
% структура параметров gNB
gNB = createNB([0, 0, HgNB], [0, 0]);
gNB.Steer = zeros(NelFull, 2);
ueRxPwr = zeros(2, N);
ueCoord = [Due, 0, 0];         % координаты UE
gNBcoords = [gNB(:).Coords].'; % массив координат gNB
distSpaceT = sqrt(sum((gNB.Coords-ueCoord.').^2)); % расстояние gNB-UE
BW = 2*atan2d(stdCoords,distSpaceT);               % ширина ДН
stAng = 0;                                         % направление максимума

% вектор, задающий направление из gNB в UE в глобальной СК x,y,z
diffCoord = ueCoord - gNBcoords;
% вектор, задающий направление из gNB в UE 
% в локальной СК АР gNB (т.е. с учетом положения АР gNB)
dirVect = gNB.AntOrient.'*diffCoord.';
% расчет углов ухода от gNB к UE
azAng = rad2deg(atan2(dirVect(2), dirVect(1)));
elAng = rad2deg(atan2(dirVect(3), sqrt(sum(dirVect(1:2).^2))));
% расчет вектора направляющих коэффициентов АР gNB
gNB.Steer(:,1) = getAntPatternSteer(antElPos, f, azAng, elAng)/NelFull;

% расчет вектора весовых коэфф. АР
[w, azAngP, antPattP] = beamshapingWeight(win_type, BW, stAng, Nel, sc);
% применять коэффициенты w если требуемая ширина ДН не меньше, 
% чем минимальная теоретическая ширина 0.891*lamb/Nel/da
if (BW*sc <= rad2deg(0.891*lamb/Nel/da) || any(isnan(w)) || sum(w) == 0)
    w = gNB.Steer(:,1);
else
    if (antType == 1)
        % расчет вектора весовых коэфф. планарной АР для вертикальных АЭ
        w = repmat(w, Nel, 1)/Nel;
        w = w(:);
    end
end
gNB.Steer(:,2) = w;

for i=1:N % цикл по числу точек расчета
    % внесение ошибки в оценку координат UE согласно stdCoords
    ueCoordErr = ueCoord;
    ueCoordErr(1:2) = ueCoordErr(1:2) + stdCoords*randn(size(ueCoord(1:2)));
    % вектор, задающий направление из gNB в UE в глобальной 
    % системе координат x,y,z с учетом ошибки координат UE
    diffCoordT = ueCoordErr - gNBcoords;
    % вектор, задающий направление из gNB в UE в системе
    % координат АР gNB (т.е. с учетом положения АР gNB)
    dirVectT = gNB.AntOrient.'*diffCoordT.';
    % расчет углов отправки от gNB к UE
    azAngT = rad2deg(atan2(dirVectT(2), dirVectT(1)));
    elAngT = rad2deg(atan2(dirVectT(3), sqrt(sum(dirVectT(1:2).^2))));

    % расчет принимаемой мощности от обслуживающей gNB 
    % с учетом диаграммообразования на gNB (без учета дальности) 
    % без/с управления шириной ДН
    gNBpwr = [getAntPatternG(antElPos, f, ...
        azAngT, elAngT, gNB.Steer(:,1), 0).^2;...
        getAntPatternG(antElPos, f, ...
        azAngT, elAngT, gNB.Steer(:,2), 0).^2];
    % расчет расстояния от UE до gNB
    diffCoord = ueCoordErr - gNBcoords;
    distSpace = sqrt(sum(diffCoord.^2,2));
    % расчет мощности принимаемой UE от gNB с учетом дальности; 
    % потери рассчитываются по модели затухания в свободном
    % пространстве для случая без/с управлением шириной ДН
    gNBpwr(isnan(gNBpwr)) = gNBpwr(1);
    gNBpwr = pow2db(gNBpwr) - fspl(distSpace,c/f);
    ueRxPwr(:, i) = gNBpwr;
end % for i=1:N % цикл по числу точек расчета

Zi = mean(ueRxPwr(2,:) - ueRxPwr(1,:));
fprintf('std=%.2f, Nel=%i, Scl=%.2f\n', s, nel, sc);
Zsc = [Zsc, Zi];
end
Zn = [Zn;Zsc];
end
Z{end+1} = Zn;
end


%%
figure;
[X,Y] = meshgrid(NellArr, sclArr);
for kk=1:size(Z,2)
    subplot(2,2,kk)
    if colorPlt == 1
        surf(X, Y, Z{kk}.', 'FaceColor', 'interp', 'EdgeColor','none');
        c1 = colorbar; c1.Label.String = '\DeltaP, дБ';
    else
        mm = ceil(max(Z{kk}(:)));
        % mn = ceil(min(Z{kk}(:)));
        mn = 0;
        me = (mm+mn)/2;
        lev = ceil([mn:5:me, me:2:mm]);
        % lev = ceil(mn:(mm-mn)/15:mm);
        [~,hh] = contour(X, Y, smoothdata(Z{kk}.', 'gaussian', 6), ...
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
if win_type == 0 
    sgtitle(sprintf('ДН Гаусса %s', antTypeCmt(antType)));
elseif win_type == 2
    sgtitle(sprintf('Прямоугольная ДН %s', antTypeCmt(antType)));
end   
fg = gcf;
if saveFig  == 1
exportgraphics(fg, strcat('pics/', picName(ww), '.tif'),'Resolution',300);
end
end

