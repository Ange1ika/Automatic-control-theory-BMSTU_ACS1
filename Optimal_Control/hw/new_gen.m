% Задаем параметры системы
sys.P = 450; 
sys.mu = 4902.8;
sys.K = 1e-3;
sys.tf = 200; 
sys.max_gen = 20; 
sys.pop_size = 10;
sys.mutation_rate = 0.3;
sys.crossover_rate = 0.8;
sys.T_theta = 10; 
sys.max_u = pi+0.01;
sys.r_moon = 1737400;

lb = [-100000000; -10000000; -2000; -200; -10; -3*pi];   
ub = [100000000; 10000000; 200; 200; 10; pi];

sys.x0 = [-2000; 1755100; 0; -1690; 1200; 0];
sys.xf = [50000; 1755050; -2; 5; 1177.5; 0.01];

function u = control_limit(lambda6, sys)
    u = max(-sys.max_u, min(sys.max_u, -lambda6 / (2 * sys.T_theta)));
end


function [J, y] = gaFunct(y, sys)
    options = odeset('Events', @(t, y) stop_conditions(t, y, sys));
    [t, y] = ode45(@(t, y) dynamics_x(t, y, sys), [0 sys.tf], y, options);
    J = my_J(y, sys);
end

function J = my_J(y, sys)

    x = y(1:6); 
    lambda = y(7:12); 
    
    u = control_limit(lambda(6), sys);
    
    norm_x = [1e4, 1e7, 1e1, 1e3, 1e3, 1]; 
    weights = [1, 1, 1, 1e5, 1, 1e3]; 
    J = 0;
    for i = 1:6
    
        J_components = sqrt(weights(i) * ((x(i) - sys.xf(i)) / norm_x(i))^2);
        J = J + J_components;
    end
end


FitnessFunction = @(lambda) gaFunct([sys.x0; lambda'], sys);

%% Структура с параметрами генетического алгоритма 
stGAopts = optimoptions("ga");
stGAopts.ConstraintTolerance = 1e-8;
stGAopts.PopulationSize = 50;
stGAopts.CrossoverFraction = 0.4;
stGAopts.EliteCount = 0.2 * stGAopts.PopulationSize;
stGAopts.FitnessLimit = -Inf;
stGAopts.FunctionTolerance = 1e-9;
stGAopts.MaxGenerations = 1000;
stGAopts.MigrationDirection = "both";
stGAopts.MaxStallGenerations = 10;
stGAopts.MaxTime = 60 * 90;
stGAopts.Display = "iter";
stGAopts.UseParallel = false;

% Поиск оптимального lambda
numOfVars = 6; 
%disp("lambda_opt");
lambda_opt = ga(FitnessFunction, numOfVars, [], [], [], [], lb, ub, [], stGAopts);

% Вывод оптимального результата
disp('Оптимальные параметры lambda:');
disp(lambda_opt);

% Условие остановки по высоте
function [value, isterminal, direction] = stop_conditions(t, y, sys)
    x = y(1:2);
    r = sqrt(x(1)^2 + x(2)^2);
    value = r - (sys.r_moon + 0.05);
    isterminal = 1;
    direction = -1;
end


options = odeset('Events', @(t, y) stop_conditions(t, y, sys));
[t, y] = ode45(@(t, y) dynamics_x(t, y, sys), [0 sys.tf], [sys.x0; lambda_opt'], options);

visualize_results(t, y, sys);

%% визуализация
function visualize_results(t, sol, sys)

    % Графики параметров
    figure('Color', 'white', 'Position', [100, 100, 1200, 600]); 
    tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    labels = {'x (Положение по оси X)', 'y (Положение по оси Y)', 'Vx (Скорость по оси X)', ...
              'Vy (Скорость по оси Y)', 'm (Масса)', 'theta (Угол)'};
    line_colors = lines(6); 

    for i = 1:6
        nexttile;
        plot(t, sol(:, i), 'LineWidth', 4, 'Color', line_colors(i, :)); 
        hold on;

        plot(t(1), sol(1, i), 'go', 'MarkerSize', 10, 'LineWidth', 2); 
        plot(t(end), sol(end, i), 'ro', 'MarkerSize', 10, 'LineWidth', 2); 
        grid on;
        title(labels{i}, 'FontSize', 14, 'FontWeight', 'bold'); 
        xlabel('Время, с', 'FontSize', 12); 
        ylabel(labels{i}, 'FontSize', 12); 
        ax = gca;
        ax.FontSize = 12; 
        ax.LineWidth = 1.5; 
        ax.GridLineStyle = '--';
    end

    % График высоты над поверхностью Луны
    r = sqrt(sol(:, 1).^2 + sol(:, 2).^2) - sys.r_moon;
    figure('Color', 'white', 'Position', [100, 100, 600, 400]);
    plot(t, r, 'LineWidth', 2, 'Color', 'b');
    hold on;
    yline(50, '--r', 'LineWidth', 1.5); % Линия ограничения
    grid on;
    title('Высота над поверхностью Луны', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Время, с', 'FontSize', 12);
    ylabel('Высота, м', 'FontSize', 12);
    ax = gca;
    ax.FontSize = 12;
    ax.LineWidth = 1.5;
    ax.GridLineStyle = '--';
    figure('Color', 'white');
    hold on;
    theta = linspace(0, 2*pi, 100);
    x_moon = sys.r_moon * cos(theta);
    y_moon = sys.r_moon * sin(theta);

    % Отрисовка Луны
    fill(x_moon, y_moon, [0.8, 0.8, 0.8], 'EdgeColor', 'none'); % Луна
    plot(sol(:, 1), sol(:, 2), 'LineWidth', 2, 'Color', 'b'); % Траектория

    plot(sol(1, 1), sol(1, 2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(sol(end, 1), sol(end, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

    grid on;
    axis equal;
    title('Траектория посадки', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('x, м', 'FontSize', 12);
    ylabel('y, м', 'FontSize', 12); 
    ax = gca;
    ax.FontSize = 12;
    ax.LineWidth = 1.5; 
    ax.GridLineStyle = '--'; 

    % График траектории
    figure('Color', 'white');
    hold on;
    theta = linspace(0, 2*pi, 100);
    x_moon = sys.r_moon * cos(theta);
    y_moon = sys.r_moon * sin(theta);

    % Отрисовка Луны
    fill(x_moon, y_moon, [0.8, 0.8, 0.8], 'EdgeColor', 'none'); % Луна
    plot(sol(:, 1), sol(:, 2), 'LineWidth', 2, 'Color', 'b'); % Траектория

    plot(sol(1, 1), sol(1, 2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(sol(end, 1), sol(end, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

    grid on;
    axis equal;
    title('Траектория посадки', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('x, м', 'FontSize', 12);
    ylabel('y, м', 'FontSize', 12); 
    ax = gca;
    ax.FontSize = 12;
    ax.LineWidth = 1.5; 
    ax.GridLineStyle = '--'; 
end


%% уравнения динамики + сопряженные переменные
function dydt = dynamics_x(t, y, sys)
    x = y(1:6); 
    lambda = y(7:12); 
    
  
    phi = atan2(x(2), x(1));
    
    %u = max(-sys.max_u, min(sys.max_u, u));
    r = sqrt(x(1)^2 + x(2)^2 + 1e-10);
    
    
    dx1 = x(3);
    dx2 = x(4);
    

    dx3 = -sys.mu/(r^2) * cos(phi) + (sys.P * sin(x(6)) * cos(phi)) / x(5) + ...
        (sys.P * cos(x(6)) * sin(phi)) / x(5);

    dx4 = -sys.mu/(r^2) * sin(phi) - (sys.P * cos(x(6)) * cos(phi)) / x(5) + ...
        (sys.P * sin(x(6)) * sin(phi)) / x(5);
    
    dx5 = -sys.K * sys.P;
    dx6 = -2 * lambda(6)/2/sys.T_theta^2 - x(6)/sys.T_theta;

    % Уравнения для сопряженных переменных
    dlambda1 = -lambda(3) * (-sys.mu * (-2 * x(1)^2 + x(2)^2) / r^5 ...
                 - (sys.P * cos(x(6)) * x(2) * x(1)) / (x(5) * r^3) ...
                 - (sys.P * sin(x(6)) * x(2)^2) / (x(5) * r^3)) ...
              - lambda(4) * (-sys.mu * x(1) * x(2) / r^5 ...
                 - (sys.P * sin(x(6)) * x(1) * x(2)) / (x(5) * r^3) ...
                 - (sys.P * cos(x(6)) * x(2)^2) / (x(5) * r^3));

    dlambda2 = -lambda(3) * (-sys.mu * x(1) * x(2) / r^5 ...
                 - (sys.P * sin(x(6)) * x(1) * x(2)) / (x(5) * r^3) ...
                 - (sys.P * cos(x(6)) * x(1)^2) / (x(5) * r^3)) ...
              - lambda(4) * (-sys.mu * (-2 * x(2)^2 + x(1)^2) / r^5 ...
                 - (sys.P * cos(x(6)) * x(2) * x(1)) / (x(5) * r^3) ...
                 - (sys.P * sin(x(6)) * x(1)^2) / (x(5) * r^3));

    dlambda3 = -2 * x(3) - lambda(1);
    dlambda4 = -2 * x(4) - lambda(2);
    dlambda5 = -(-lambda(3) * (-sys.P * cos(x(6)) * x(2) + sys.P * sin(x(6)) * x(1)) / ...
                 (x(5)^2 * r) ...
                 + lambda(4) * (-sys.P * sin(x(6)) * x(2) + sys.P * cos(x(6)) * x(1)) / ...
                 (x(5)^2 * r));
    dlambda6 = -(lambda(3) * ((sys.P * cos(x(6)) * x(2) - sys.P * sin(x(6)) * x(1)) / (x(5) * sqrt(x(1)^2 + x(2)^2)))) - ...
        lambda(4) * ((-sys.P * sin(x(6)) * x(1) + sys.P * cos(x(6)) * x(2)) / (x(5) * r)) + lambda(6) / sys.T_theta; 
    dydt = [dx1; dx2; dx3; dx4; dx5; dx6; dlambda1; dlambda2; dlambda3; dlambda4; dlambda5; dlambda6];
end

