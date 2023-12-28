c = physconst('LightSpeed');
anim = 0;       % 1 - вкл. анимацию работы
f = 30e9;       % несущая в диапазоне ММВ, Гц
lamb = c/f;     % длина волны, м
da = 0.5*c/f;   % расстояние между элементами АР
snrThr = 10;    % пороговое отношения сигнал/помеха для отображения карты
useAntUE = 0;   % использовать ДН АР на UE (только режим antPattCntrl=0)

backLobe = 1; % подавление обратного лепестка ДН (логично только для URA)

% отображение ДН в 2D при выставленном флаге анимации (anim = 1)
antPlot2D = 1; % (рекомендуется выставлять в 1)
% Шаг отрисовки ДН при anim = 1
if (antPlot2D == 0)
    angStep = 5;
else
    angStep = 1;
end

% период процедуры диаграммообразования, с;
Ta = 0; % при Ta = 0 ДН формируется на каждом шаге расчета


antElPos = createAnt(antType, Nel, da); % формирование АР
NelFull = size(antElPos, 1);            % полное число АЭ

% выбор сценария:
% 1 - gNB расположены на одной стороне относительно траектории движения UE
% 2 - gNB расположены по разные стороны относительно траектории движения UE
[gNB, ueNode, d, T, v] = createScenarion(1);
Nd = length(d);

Nnb = length(gNB);    % число gNB
Nue = length(ueNode); % число UE
% число точек расчета (число точек координат траектории UE)
N = length(ueNode(1).Trajectory(:,1));
trajArray = [ueNode.Trajectory];       % массив координат UE [N x 3*Nue]
gNBcoords = [gNB(:).Coords].';         % массив координат gNB [Nnb x 3]


% инициализация массивов для хранения углов отправки от каждой gNB 
% к каждой UE (нужно для устранении ошибки чтения несуществующего 
% массива при некоторых настройках модели)
azAng = zeros(Nue, Nnb);
elAng = zeros(Nue, Nnb);
azAngUE = zeros(Nue, Nnb);
elAngUE = zeros(Nue, Nnb);

%% ЦИКЛ ПО ЧИСЛУ ТОЧЕК РАСЧЕТА
for i=1:N % цикл по числу точек расчета
    % массив координат всех UE для i-й точки расчета
    ueCoordsi = reshape(trajArray(i, :).', 3, Nue).';
    % внесение ошибки в оценку координат UE согласно stdCoords
    ueCoordsiErr = ueCoordsi;
    ueCoordsiErr(:,1:2) = ueCoordsiErr(:,1:2) + stdCoords*randn(2,2);
    if (mod(i, round(Ta/T)) == 1 || Ta == 0)
        % инициализация массивов для хранения углов ухода 
        % от каждой gNB к каждой UE
        azAng = zeros(Nue, Nnb);  % угол ухода по азимуту
        elAng = zeros(Nue, Nnb);  % угол ухода по углу места
        azAngT = zeros(Nue, Nnb); % истинный угол ухода по азимуту
        elAngT = zeros(Nue, Nnb); % истинный угол ухода по углу места
        % инициализация массивов для хранения углов прихода 
        % для каждой UE от каждой gNB (используется при наличии АР на UE)
        azAngUE = zeros(Nue, Nnb);
        elAngUE = zeros(Nue, Nnb);
        for j=1:Nue % цикл по числу UE
            % вектор, задающий направление из gNB в UE в глобальной СКx,y,z
            diffCoord = ueCoordsiErr(j,:) - gNBcoords;
            diffCoordT = ueCoordsi(j,:) - gNBcoords;
            for n=1:Nnb % цикл по числу gNB
                % вектор, задающий направление из gNB в UE в локальной 
                % системе координат АР gNB (т.е. с учетом положения АР gNB)
                dirVect = gNB(n).AntOrient.'*diffCoord(n,:).';
                dirVectT = gNB(n).AntOrient.'*diffCoordT(n,:).';
                % расчет углов ухода от n-й gNB к j-й UE
                azAng(j, n) = rad2deg(atan2(dirVect(2), dirVect(1)));
                elAng(j, n) = rad2deg(atan2(dirVect(3), ...
                    sqrt(sum(dirVect(1:2).^2))));
                azAngT(j, n) = rad2deg(atan2(dirVectT(2), dirVectT(1)));
                elAngT(j, n) = rad2deg(atan2(dirVectT(3), ...
                    sqrt(sum(dirVectT(1:2).^2))));
                % расчет углов ухода для АР UE
                if ( useAntUE == 1)
                    % вектор, задающий направление из UE в gNB в локальной
                    % системе координат АР UE,т.е. с учетом положения АР UE
                    dirVect = -ueNode(j).AntOrient.'*diffCoord(n,:).';
                    % расчет углов ухода от j-й UE к n-й gNB
                    azAngUE(j, n) = rad2deg(atan2(dirVect(2), dirVect(1)));
                    elAngUE(j, n) = rad2deg(atan2(dirVect(3), ...
                        sqrt(sum(dirVect(1:2).^2))));
                end
            end
        end
        % для адаптивной схемы, т.к. она не использует информацию 
        % о местоположении UE
        if (antPattCntrl == 2)
            azAng = azAngT;
            elAng = elAngT;
        end
        
        for n=1:Nnb % цикл по числу gNB
            cVect = zeros(NelFull, Nue);
            gVect = zeros(1, Nue);
            for j=1:Nue % цикл по числу UE       
                % расчет направляющего вектора 
                % фазового распределения АР n-й gNB к j-й UE
                cVect(:,j) = ...
                 getAntPatternSteer(antElPos, f, azAng(j, n), elAng(j, n));
                
                % сохранения данных о обслуживаемой UE (номер обслуживающей 
                % gNB указывается в параметре servgNB для каждой UE)
                if (n == ueNode(j).servgNB)
                    % расчет вектора направляющих коэффициентов
                    % АР n-й gNB обслуживающей j-ю UE 
                    gNB(n).Steer = getAntPatternSteer(antElPos, f, ...
                        azAng(j, n), elAng(j, n))/NelFull;
                    % запись 1 в векторе g (см. формулу (10)) по индексу
                    % j-й UE для выставления максимума в этом направлении
                    gVect(j) = 1;
                    % сохранение информации о местоположении UE в полярных
                    % координатах (азимут, угол места, расстояние)
                    gNB(n).UEPolarCoord = [azAng(j, n), elAng(j, n),...
                    sqrt(sum((ueCoordsiErr(j,:).' - gNB(n).Coords).^2))];
                end

                % расчет вектора направляющих коэффициентов
                % АР j-й gNB работающей с n-й gNB
                if (n == ueNode(j).servgNB && useAntUE == 1)
                    ueNode(j).Steer = getAntPatternSteer(antElPos, f, ...
                        azAngUE(j,n), elAngUE(j,n));
                end
            end
            % расчет вектора направляющих коэфф. АР n-й gNB обслуживающей
            % j-ю UE при схеме управления ДН отличной от antPattCntrl=0
            switch antPattCntrl
                % LCMV
                case 1
                    % расчет вектора коэффициентов по алгоритму LCMV для
                    % управления максимумом и нулем ДН (см. формулу (13))
                    gNB(n).Steer = cVect*pinv(cVect'*cVect)*gVect';
                    % нормировка коэффициентов
                    gNB(n).Steer = ...
                        gNB(n).Steer/max(abs(gNB(n).Steer))/NelFull;
                % SMI
                case 2
                    K = 100; % длина преамбулы
                    % направляющий вектор АР для полезного сигнала
                    sv_s = cVect(:,[ueNode.servgNB] == n);
                    % направляющий вектор АР для помехи
                    sv_j = cVect(:,~([ueNode.servgNB] == n));
                    % симуляция приема сигнала и помехи с заданных
                    % направлений
                    rx_sum = ones(K,1)*sv_s.' + ...
                        (randn(K,1) + 1i*randn(K,1))/sqrt(2)*sv_j.';
                    % корреляционная матрица сигналов АЭ (см. формулу (22))
                    Rxx = rx_sum'*rx_sum;
                    % корреляционный вектор сигналов АЭ 
                    % с опорной преамбулой (см. формулу (22))
                    rxd = rx_sum'*ones(K,1);
                    % расчет вектора коэфф. по алгоритму SMI 
                    % для подавления помехи (см. формулу (22))
                    gNB(n).Steer = conj(pinv(Rxx)*rxd);
                % beamshaping
                case 3
                    % ширина ДН, рассчитанная по СКО оценки координат UE
                    % (см. формулу (41))
                    BW = 2*atan2d(stdCoords,gNB(n).UEPolarCoord(3));
                    % расчет вектора весовых коэфф. АР
                    w = beamshapingWeight(win_type, BW, ...
                        gNB(n).UEPolarCoord(1), Nel, sc);
                    % применять коэфф. w если требуемая ширина ДН:
                    % не меньше, чем минимальная теоретическая ширина
                    % равная 0.891*lamb/Nel/da; амплитуда суммы коэфф.>0.1
                    if ~(BW*sc <= rad2deg(0.891*lamb/Nel/da) ||...
                            any(isnan(w)) || sum(w) < 0.1)
                        if (antType == 1)
                            % расчет вектора весовых коэффициентов
                            % планарной АР для вертикальных АЭ
                            antElPosElev = [zeros(Nel,1), ...
                                antElPos(1:Nel,2), zeros(Nel,1)];
                            wEl = getAntPatternSteer(antElPosElev, ...
                                f, gNB(n).UEPolarCoord(2), 0)/Nel;
                            % вектор весовых коэффициентов планарной АР,
                            % полученный перемножением горизонтальных и 
                            % вертикальных коэффициентов
                            w = w*wEl.'/Nel;
                            w = w(:);
                        end                        
                        gNB(n).Steer = w;
                    end
            end % switch antPattCntrl
        end % for n=1:Nnb % цикл по числу gNB       
    end % if (mod(i, round(Ta/T)) == 1 || Ta == 0)
    
    % массив для временного хранения значений 
    % принимаемой мощности на UE от каждой gNB
    gNBpwr = zeros(Nnb, 1);

    % расчет отношения принимаемой мощности от обслуживающей gNB 
    % к мощности, принимаемой от соседней gNB для каждой UE
    for j=1:Nue % цикл по числу UE
        for n=1:Nnb % цикл по числу gNB
            % расчет коэфф усиления приемной антенны UE
            if ( useAntUE == 1 && j == 1)
                % КУ при наличии АР на UE, луч которой 
                % направлен на обслуживающую gNB
                gUE = getAntPatternG(antElPos, f, ...
                    azAngUE(n), elAngUE(n), ueNode(1).Steer, backLobe).^2;
            else
                % КУ без АР на UE
                gUE = 1;
            end
            % расчет мощности принимаемой на j-й UE от n-й gNB с учетом
            % диаграммообразования на UE и gNB (без учета дальности)
            gNBpwr(n) = gUE*getAntPatternG(antElPos, f, ...
                azAngT(j, n), elAngT(j, n), gNB(n).Steer, backLobe).^2;
        end
        % расчет расстояния от j-й UE до каждой gNB
        diffCoord = ueCoordsi(j,:) - gNBcoords;
        distSpace = sqrt(sum(diffCoord.^2,2));
        % расчет мощности принимаемой на j-й UE от каждой gNB с учетом
        % дальности (потери рассчитываются по модели затухания в свободном
        % пространстве)
        gNBpwr = pow2db(gNBpwr) - fspl(distSpace,c/f);
        % расчет отношения принимаемой мощности от обслуживающей gNB 
        % к мощности принимаемой от соседней gNB для j-й UE
        ueNode(j).SNR(i) = gNBpwr(ueNode(j).servgNB) - ...
            sum(gNBpwr(1:end ~= ueNode(j).servgNB));
    end
    
    % обновление отображения положения UE и ДН АР gNB/UE 
    % для текущего шага расчета
    if (anim == 1)
        for ip=1:Nnb
            % обновление ДН gNB
            [x, y, z] = antPattPlot(antElPos, f, ...
                gNB(ip), angStep, antPattScl, backLobe, antPlot2D);
            gNBptrnPlot(ip).XData = x + gNB(ip).Coords(1);
            gNBptrnPlot(ip).YData = y + gNB(ip).Coords(2);
            if (antPlot2D == 0)
                gNBptrnPlot(ip).ZData = z + gNB(ip).Coords(3);
            end
            % обновление ДН UE
            if (useAntUE == 1)
                [x, y, z] = antPattPlot(antElPos, f, ...
                    ueNode(ip),angStep, antPattSclUe, backLobe, antPlot2D);
                ueptrnPlot(ip).XData = x + ueCoordsi(ip,1);
                ueptrnPlot(ip).YData = y + ueCoordsi(ip,2);
                if (antPlot2D == 0)
                    ueptrnPlot(ip).ZData = z + ueCoordsi(ip,3);
                end
            end
            % обновление вектора направления от gNB на UE
            ueDirPlot(ip).XData = [gNB(ip).Coords(1); ueCoordsi(ip,1)];
            ueDirPlot(ip).YData = [gNB(ip).Coords(2); ueCoordsi(ip,2)];
            ueDirPlot(ip).ZData = [gNB(ip).Coords(3); ueCoordsi(ip,3)];
            % обновление вектора направления от gNB на  соседнюю UE
            ueDirPlot2(ip).XData = [gNB(ip).Coords(1); ...
                ueCoordsi(indSnoi(ip),1)];
            ueDirPlot2(ip).YData = [gNB(ip).Coords(2); ...
                ueCoordsi(indSnoi(ip),2)];
            ueDirPlot2(ip).ZData = [gNB(ip).Coords(3); ...
                ueCoordsi(indSnoi(ip),3)];
            % обновление положения UE
            uePlot(ip).XData = ueCoordsi(ip,1);
            uePlot(ip).YData = ueCoordsi(ip,2);
            uePlot(ip).ZData = ueCoordsi(ip,3);
            ueText(ip).Position = [ueCoordsi(ip,1)+2, ...
                                   ueCoordsi(ip,2), ...
                                   ueCoordsi(ip,3)+5];
        end  % for ip=1:Nnb    
        pause(0.05)
    end % if (anim == 1)
end % for i=1:N % цикл по числу точек расчета

%%
% подготовка массивов координат и массива значений 
% отношения сигнал/помеха 1-й UE для отображения
X = reshape(ueNode(1).Trajectory(:,1), [], Nd).';
Y = reshape(ueNode(1).Trajectory(:,2), [], Nd).';
Z = reshape(ueNode(1).SNR, [], Nd).';
