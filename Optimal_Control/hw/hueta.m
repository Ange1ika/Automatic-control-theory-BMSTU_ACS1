% Очистка рабочего пространства
clc; clear; close all;

% Задание параметров системы
sys.P = 450; 
sys.mu = 4902.8;
sys.K = 0.01; 
sys.T_theta = 10;

% Начальные условия состояния
x1_0 = 0;         % Начальная x позиция
x2_0 = 18000;     % Начальная y позиция
x3_0 = 0;         % Начальная скорость по x
x4_0 = 1690;      % Начальная скорость по y
x5_0 = 1200;      % Начальная масса
x6_0 = 0;         % Начальный курс

sys.x0 = [x1_0; x2_0; x3_0; x4_0; x5_0; x6_0];

% Конечные условия состояния
sys.xf = [50; 1737050; -2; 5; 975; 10];

% Начальное приближение для сопряжённых переменных и конечного времени
lambda0_guess = zeros(6, 1);
tf_guess = 100; % Предположительное терминальное время

% Вектор переменных оптимизации
decision_vars_guess = [lambda0_guess; tf_guess];

% Границы для сопряжённых переменных и терминального времени
lb = [-1e3 * ones(6, 1); 500];   % Минимум: tf >= 10 секунд
ub = [1e3 * ones(6, 1); 500];   % Максимум: tf <= 500 секунд

% Определение целевой функции для генетического алгоритма
cost_function = @(decision_vars) compute_cost(decision_vars, sys);

% Настройки генетического алгоритма
options = optimoptions('ga', 'Display', 'iter', 'PopulationSize', 100, ...
    'MaxGenerations', 50, 'FunctionTolerance', 1e-6);

% Оптимизация с использованием GA
[decision_vars_opt, fval] = ga(cost_function, 7, [], [], [], [], lb, ub, [], options);

% Извлечение оптимальных параметров
lambda0_opt = decision_vars_opt(1:6);
tf_opt = decision_vars_opt(7);

% Решение системы ODE с оптимальными параметрами
tspan = [0, tf_opt];
y0 = [sys.x0; lambda0_opt(:)];
[t, y] = ode45(@(t, y) dynamics(t, y, sys), tspan, y0);

% Выделение переменных состояния и сопряжённых переменных
x = y(:, 1:6);
lambda = y(:, 7:12);

% Расчёт управления, лагранжиана и гамильтониана
u = zeros(length(t), 1);
L = zeros(length(t), 1);
H = zeros(length(t), 1);
for i = 1:length(t)
    % Оптимальное управление
    u(i) = compute_control(y(i, 1:6), y(i, 7:12), sys.T_theta);
    
    % Лагранжиан
    L(i) = sum(y(i, 1:6).^2) + u(i)^2; 
    
    % Гамильтониан
    H(i) = compute_hamiltonian(y(i, 1:6), y(i, 7:12), u(i), sys);
end

% Интегральная стоимость
J = trapz(t, L);

% Визуализация результатов
visualize_results(t, x, lambda, u, H);

% --- Функции ---

function visualize_results(t, x, lambda, u, H)
    % Визуализация состояния
    figure;
    tiledlayout(6, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:6
        nexttile;
        plot(t, x(:, i), 'LineWidth', 2);
        grid on; 
        title(sprintf('x%d vs time', i));
    end

    % Управление
    figure;
    plot(t, u, 'LineWidth', 2);
    grid on;
    title('Control (u) vs time');
    xlabel('Time');
    ylabel('u');

    % Гамильтониан
    figure;
    plot(t, H, 'LineWidth', 2);
    grid on;
    title('Hamiltonian (H) vs time');
    xlabel('Time');
    ylabel('H');
end

function cost = compute_cost(decision_vars, sys)
    lambda0 = decision_vars(1:6);
    tf = decision_vars(7);

    % Интеграция системы с текущими параметрами
    tspan = [0 tf];
    y0 = [sys.x0; lambda0(:)];
    [~, y] = ode45(@(t, y) dynamics(t, y, sys), tspan, y0);

    % Конечное состояние
    final_state = y(end, 1:6)';

    % Стоимость — норма отклонения от терминальных условий
    cost = norm(final_state - sys.xf)^2;
end

function dydt = dynamics(t, y, sys)
    % Переменные состояния и сопряжённые переменные
    x = y(1:6); 
    lambda = y(7:12);

    % Дополнительные параметры
    phi = atan2(x(2), x(1));
    u = compute_control(x, lambda, sys.T_theta);
    r = sqrt(x(1)^2 + x(2)^2 + 1e-10);

    % Дифференциальные уравнения
    dx = zeros(6, 1);
    dx(1) = x(3);
    dx(2) = x(4);
    dx(3) = -sys.mu / r^2 * sin(phi) - (sys.P * cos(x(6)) * cos(phi)) / x(5) + ...
            (sys.P * sin(x(6)) * sin(phi)) / x(5);
    dx(4) = -sys.mu / r^2 * cos(phi) + (sys.P * sin(x(6)) * cos(phi)) / x(5) + ...
            (sys.P * cos(x(6)) * sin(phi)) / x(5);
    dx(5) = -sys.K * sys.P;
    dx(6) = (u - x(6)) / sys.T_theta;

    % Уравнения сопряжённых переменных
    dlambda = zeros(6, 1);
    dlambda(1) = -lambda(3) * (sys.mu * x(1) / r^3) + ...
                 lambda(4) * (sys.mu * x(2) / r^3);
    dlambda(2) = -lambda(3) * (sys.mu * x(2) / r^3) - ...
                 lambda(4) * (sys.mu * x(1) / r^3);
    dlambda(3) = -lambda(1);
    dlambda(4) = -lambda(2);
    dlambda(5) = lambda(3) * (sys.P * cos(x(6)) / x(5)^2) - ...
                 lambda(4) * (sys.P * sin(x(6)) / x(5)^2);
    dlambda(6) = lambda(3) * (sys.P * sin(x(6)) / x(5)) - ...
                 lambda(4) * (sys.P * cos(x(6)) / x(5));

    % Собираем уравнения
    dydt = [dx; dlambda];
end

function u = compute_control(x, lambda, T_theta)
    u = -lambda(6) / (2 * T_theta);
end

function H = compute_hamiltonian(x, lambda, u, sys)
    dx = dynamics(0, [x; lambda], sys);
    H = dot(lambda, dx(1:6)) + u^2;
end
