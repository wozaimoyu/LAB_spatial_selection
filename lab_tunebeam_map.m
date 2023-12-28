% Построение рис. 8, 9

clear all; close all; clc;
mkdir pics

% Цветовая гамма рисунков
% 1 - цветные рисунки, 0 - ч\б рисунки
colorPlt = 0;
% автоматическое сохранение рисунков
saveFig = 1;



stdCoordsArr = 10;  % СКО оценки координат UE по осям x, y, z в метрах
Nel = 20;                               % число АЭ в одном измерении
% выбор типа антенной решетки:
% 1 - планарная АР
% 2 - линейная АР
% 3 - круговая АР (не поддерживается)
antType = 1;


picName = ["рис8", "рис9"];


% выбор формы ДН для алгоритма управления шириной ДН (antPattCntrl = 3)
% 0 - окно Гаусса
% 1 - окно приподнятого косинуса
% 2 - прямоугольное окно
win_typeArr = [2, 0];
antTypeCmt = ["Планарная АР", "Линейная АР", "Круговая АР"];

ueRxPwrPlt = [];
figNumber = 1;
for s=stdCoordsArr
subFigNumber = 1;
for ww=win_typeArr
c = physconst('LightSpeed');
f = 30e9;       % несущая в диапазоне ММВ, Гц
lamb = c/f;     % длина волны, м
da = 0.5*c/f;   % расстояние между элементами АР
snrThr = 10;    % пороговое отношения сигнал/помеха для отображения карты
stdCoords = s; % СКО оценки координат x, y, z


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

scl = 1.5; % множитель ширины ДН (см. формулу (42), (43))
% расчет вектора весовых коэфф. АР
[w, azAngP, antPattP] = beamshapingWeight(ww, BW, stAng, Nel, scl);
% применять коэффициенты w если требуемая ширина ДН не меньше, 
% чем минимальная теоретическая ширина 0.891*lamb/Nel/da
if (BW*scl <= rad2deg(0.891*lamb/Nel/da) || any(isnan(w)))
    w = gNB.Steer(:,1);
else
    if (antType == 1)
        % расчет вектора весовых коэфф. планарной АР для вертикальных АЭ
        w = repmat(w, Nel, 1)/Nel;
        w = w(:);
    end
end
gNB.Steer(:,2) = w;
rng('default')
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

% рассчет ДН для случая без/с управлением шириной ДН
azAngPatt = -90:0.1:89;
g = zeros(1, length(azAngPatt));
gDef = zeros(1, length(azAngPatt));
for i=1:length(azAngPatt)
    g(i) = getAntPatternG(antElPos, f, azAngPatt(i), 0, w, 0);
    gDef(i) = getAntPatternG(antElPos, f, ...
        azAngPatt(i), 0, ones(size(w))/NelFull, 0);
end
gNorm = g/max(g);
% рассчет ширины ДН по уровню -3 дБ
[~,ind] = min(abs(gNorm - 1/sqrt(2)));
hpbw = 2*abs(azAngPatt(ind));

% расчет координат кривых ДН для отображения
% повернуто на 90 градусов для лучшей визуализации
alphPatt = azAngPatt + 90;
xPatt = cosd(alphPatt).*gNorm*(distSpaceT);
yPatt = sind(alphPatt).*gNorm*(distSpaceT);
xPattDef = cosd(alphPatt).*gDef*(distSpaceT);
yPattDef = sind(alphPatt).*gDef*(distSpaceT);

% отображение ДН на карте вероятности положения UE
figure(figNumber); 
subplot(1,2,subFigNumber); 
if colorPlt == 1
    % построение цветовой карты вероятности положения UE
    [X,Y] = meshgrid(-25:25, 0:Due+15);
    p = mvnpdf([X(:) Y(:)], ueCoord([2,1]), diag([stdCoords,stdCoords].^2));
    p = reshape(p,size(X));
    pcolor(X,Y,p);
    shading interp;
    hold on;
    plot(xPatt, yPatt, 'Color', '#D95319', 'LineWidth', 1.3);
    plot(xPattDef, yPattDef, 'Color', '#7E2F8E', 'LineWidth', 1.3);
    c = colorbar; c.Label.String = 'p(x,y)';
    axis tight;
else
    plot(xPatt, yPatt, 'k--', 'LineWidth', 1.1);
    hold on;
    plot(xPattDef, yPattDef, 'k-', 'LineWidth', 1.1);
    ax = gca;
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    ax.GridColor = [0 0 0];
    prbA = [0.1, 0.5, 0.9];
    [~, ~, ppR] = get_prob(ueCoord(2), ueCoord(1), stdCoords, prbA);
    angA = 0:360;
    for i=1:length(prbA)
        plot(ppR(i)*cosd(angA) + ueCoord(2), ...
             ppR(i)*sind(angA) + ueCoord(1), 'k-');
        text(ueCoord(2), ppR(i) + ueCoord(1), num2str(prbA(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    ylim([0, ppR(end) + ueCoord(1)+3])
end
text(ueCoord(2), ueCoord(1), 'UE', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom')
grid on; axis equal; xlabel('x, м'); ylabel('y, м');
if ww == 0 
    title('ДН Гаусса');
elseif ww == 2
    title('Прямоугольная ДН');
end
subFigNumber = subFigNumber + 1;
ueRxPwrPlt = [ueRxPwrPlt; ueRxPwr];
end
if saveFig  == 1
exportgraphics(figure(figNumber), ...
    strcat('pics/', picName(figNumber), '.tif'),'Resolution',300);
end
sgtitle([sprintf('СКО оценки координат UE %i м', s), antTypeCmt(antType)]);
figNumber = figNumber + 1;
end

figure;
if colorPlt == 1
    plot(ueRxPwrPlt([1, 4, 2],:).')
else
    plot(ueRxPwrPlt(1,:).', 'k-');
    hold on
    plot(ueRxPwrPlt(2,:).', 'kx-');
    plot(ueRxPwrPlt(4,:).', 'ko-');
    ax = gca;
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    ax.GridColor = [0 0 0];
end
grid on; xlabel('Номер точки расчета'); ylabel('P, дБ');
legend('без упр ДН', 'Прямоуг', 'Гаусс', 'Location', 'southeast');
if saveFig  == 1
exportgraphics(ax, ...
    strcat('pics/', picName(figNumber), '.tif'),'Resolution',300);
end
