x = -10:0.1:10;
y = x;

function Z = my_J(the)
    x = the(1);
    y = the(2);
    % R = sqrt(X.^2 + Y.^2) + eps;
    % Z = sin(R)./R;
    % J = - Z;

    Z =  3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ... 
   - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ... 
   - 1/3*exp(-(x+1).^2 - y.^2) 
end

function Z = my_func(x, y)

    Z =  3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ... 
   - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ... 
   - 1/3*exp(-(x+1).^2 - y.^2) 

end


[X, Y] = meshgrid(x, y);
Z = my_func(X, Y);

% GA

%% Структура с параметрами генетического алгоритма 
stGAopts = optimoptions("ga");
% Допустимая ошибка выполнения ограничений
stGAopts.ConstraintTolerance = 1e-8;
% Размер популяции, т.е. количество особей в каждом поколении
stGAopts.PopulationSize = 500;
% Доля популяции, участвующая в кроссовере (скрещивании)
stGAopts.CrossoverFraction = 0.2;
% Количество элитных особей, которые будут автоматически перенесены в следующее поколение
stGAopts.EliteCount = 0.5 * stGAopts.PopulationSize;
% Ограничение на значение целевой функции (функция приспособленности), при котором алгоритм останавливается
stGAopts.FitnessLimit = -Inf;
stGAopts.FunctionTolerance = 1e-8;
% Максимальное количество поколений, которое будет создано
stGAopts.MaxGenerations = 10000;
% Направление миграции особей между субпопуляциями (в данном случае - в обоих направлениях)
stGAopts.MigrationDirection = "both"
% Максимальное количество поколений, в течение которых значение целевой функции может оставаться неизменным (стагнация)
stGAopts.MaxStallGenerations = 5;
stGAopts.MaxTime = 60*90;
stGAopts.Display = "iter";
stGAopts.UseParallel = false; % true;

%%
numOfVars = 2;
% GA зависит только от одного вектора, но для доп структур
% анонимная функция
FitnessFunction = @(the) my_J(the);
THE = ga(FitnessFunction, numOfVars);
% the = [x, y]

% График 
h1 = figure("Units", "normalized", "OuterPosition", ...
    [0.05 0.05 0.9 0.88]);
clf;


surf(X, Y, Z);
hold on;
plot3(THE(1), THE(2), my_func(THE(1), THE(2)), "*r", "MarkerSize", 15, "LineWidth", 5);