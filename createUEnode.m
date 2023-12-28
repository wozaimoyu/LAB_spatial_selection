% функция создания структуры UE
% входные параметры:
% xPnts - начальная и конечные значения траектории UE по X, м
% yPnts - начальная и конечные значения траектории UE по Y, м
% zPnts - высота UE (постоянна на всей траектории), м
% servgNB - номер обслуживающей gNB (индекс в массиве структур gNB)
% v - скорость движения UE, м/с
% T - период измерения координат UE, с
% antDir - ориентация АР UE (направление оси симметрии) 
%          [азимут, наклон], градусы
function ue = createUEnode(xPnts, yPnts, zPnts, servgNB, v, T, antDir)
ue.servgNB = servgNB; % номер обслуживающей eNB    
% создание массива координат траектории движения UE
trajPnts = [xPnts, yPnts];
trajPnts(:,3) = zPnts;
ue.Trajectory = getTrajectory(trajPnts, v, T);    
% инициализация массива значений отношения сигнал/помеха,
% рассчитывается в каждой точке траектории UE
ue.SNR = zeros(length(ue.Trajectory(:,1)), 1);    
% матрица поворота согласно значениям antDir, используется для пересчета 
% векторов направлений из глобальных координат в систему координат АР UE
ue.AntOrient = rotz(antDir(1))*roty(-antDir(2));
% вектор направляющих коэфф. АР (поумолчанию АР не направлена)
ue.Steer = 1;
end