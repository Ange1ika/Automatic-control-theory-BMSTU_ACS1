% Метод роя частиц
%% Структура с параметрами 
stPSopts = optimoptions("particleswarm");

stPSopts.FunctionTolerance = 1e-8;
stPSopts.MaxIterations = 10000;
stPSopts.MaxTime = 60*90;
stPSopts.SwarmSize = 3; % Число частиц в рое
stPSopts.SelfAdjustmentWeight = 1.49;
stPSopts.SocialAdjustmentWeight = 1.49;
stPSopts.Display = "iter";
stPSopts.UseParallel = false; % true;
%%
% Поиск оптимального lambda
numOfVars = 2; % Размерность вектора lambda
lb = [-500,  -500];
ub = [500, 500];

x = -500:10:500;
y = x;
[X, Y] = meshgrid(x, y);
Z = my_func(X, Y);

% МРЧ
%FitnessFunction = @(the) my_J_
THE = particleswarm(@my_J, numOfVars, lb, ub, stPSopts)


h3 = figure("Units", "normalized", "OuterPosition", [0.05 0.05 0.9 0.88])
clf;

surf(X, Y, Z)
hold on;
plot3(THE(1), THE(2), my_func(THE(1), THE(2)), "rx", "LineWidth", 6, "MarkerSize", 50)

%%
function J = my_J(the)
    X = the(1);
    Y = the(2);

    J = -X.*sin(abs(X).^0.5)-Y.*sin(abs(Y).^0.5);

end


function f = my_func(X, Y)

    f = -X.*sin(abs(X).^0.5)-Y.*sin(abs(Y).^0.5);
end