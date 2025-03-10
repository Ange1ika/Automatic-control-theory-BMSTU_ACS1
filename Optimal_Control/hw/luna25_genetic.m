function luna25_genetic
    % Параметры системы
    sys.P = 450; 
    sys.mu = 4902.8;
    sys.K = 1e-3;
    sys.tf = 200; 
    sys.max_gen = 5; 
    sys.pop_size = 10;
    sys.mutation_rate = 0.3;
    sys.crossover_rate = 0.8;
    sys.T_theta = 10; 
    sys.max_u = 100 * pi;
    sys.r_moon = 1737400;
    sys.epsilon = 1e-4;
    sys.min_height = 0.05; 
    sys.norm_x = [1e7, 1e7, 1e1, 1e3, 1e3, 1]; 

    sys.x0 = [-200000; 1755100; 0; -1690; 1200; 0];
    sys.xf = [5000000; 1755050; -2; 5; 1177.5; 0.01];

    % Весовые коэффициенты
    sys.weights = [1, 1, 0.1, 1e3, 1, 0.01];

    % Границы параметров
    lb = [-100000000; -10000000; -2000; -200; -10; -3*pi];   
    ub = [100000000; 10000000; 200; 200; 10; pi];
    

    % Генерация начальной популяции
    population = initialize_population(sys.pop_size, lb, ub);

    % ГА
    for gen = 1:sys.max_gen
        
        fitness = evaluate_population(population, sys);

        % Отбор
        selected = tournament_selection(population, fitness, sys.pop_size);

        % Кроссовер
        offspring = crossover(selected, sys.crossover_rate, lb, ub);

        % Мутация
        offspring = mutate(offspring, sys.mutation_rate, lb, ub);

        % Обновление популяции
        population = offspring;

        [min_fitness, best_idx] = min(fitness);
        disp(['Поколение ', num2str(gen), ': Лучшая ошибка = ', num2str(min_fitness)]);

        if min_fitness < sys.epsilon
            break;
        end
    end

    best_lambda = population(best_idx, :);
    disp('Оптимальное значение λ:');
    disp(best_lambda);

    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, ...
                     'Events', @(t, y) height_event(t, y, sys));
    [t, sol] = ode45(@(t, y) dynamics(t, y, sys), [0 sys.tf], [sys.x0; best_lambda(:)], options);

    visualize_results(t, sol, sys);
end

%% инициализация нач популяции
function population = initialize_population(pop_size, lb, ub)

    num_params = length(lb);
   
    population = zeros(pop_size, num_params);

    for i = 1:num_params
        population(:, i) = lb(i) + (ub(i) - lb(i)) * rand(pop_size, 1);
    end
end

%% функционал качества
function fitness = evaluate_population(population, sys)
    num_individuals = size(population, 1);
    fitness = zeros(num_individuals, 1);
    
    for i = 1:num_individuals
        lambda0 = population(i, :)';
        [t, xf, valid] = simulate_dynamics(sys, lambda0);
        error = abs((xf - sys.xf) ./ sys.norm_x) .* sys.weights(:);
        fitness(i) = sum(error);
        
    end
end

%% решение ДУ

function [t, xf, valid] = simulate_dynamics(sys, lambda0)
    y0 = [sys.x0; lambda0(:)];
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, ...
                     'Events', @(t, y) height_event(t, y, sys));
    [t, y, te, ye, ie] = ode45(@(t, y) dynamics(t, y, sys), [0 sys.tf], y0, options);
    xf = y(end, 1:6)'; 
    valid = isempty(ie); % Проверка, достигнута ли высота
end

%% проверка на достижение высоты
function [value, isterminal, direction] = height_event(t, y, sys)
    r = sqrt(y(1)^2 + y(2)^2) - sys.r_moon;
    value = r - sys.min_height; 
    isterminal = 1; % Останавливаем интегратор
    direction = -1; % Проверяем только пересечение сверху вниз
end

%% отбор элитных популяций
function selected = tournament_selection(population, fitness, pop_size)
    num_individuals = size(population, 1);
    selected = zeros(pop_size, size(population, 2));
    for i = 1:pop_size
        idx = randi(num_individuals, 2, 1); % Выбор двух случайных индивидов
        if fitness(idx(1)) < fitness(idx(2))
            selected(i, :) = population(idx(1), :);
        else
            selected(i, :) = population(idx(2), :);
        end
    end
end

%% кроссинговер
function offspring = crossover(population, rate, lb, ub)
    num_individuals = size(population, 1);
    num_vars = size(population, 2);
    offspring = population;

    for i = 1:2:num_individuals-1
        if rand < rate
            % Создаем маску для выбора генов
            mask = rand(1, num_vars) > 0.5;
            parent1 = population(i, :);
            parent2 = population(i+1, :);

            % Применяем кроссовер
            offspring(i, :) = parent1 .* mask + parent2 .* ~mask;
            offspring(i+1, :) = parent2 .* mask + parent1 .* ~mask;

            % Ограничиваем значения в пределах границ для каждого элемента
            offspring(i, :) = max(lb', min(ub', offspring(i, :)));
            offspring(i+1, :) = max(lb', min(ub', offspring(i+1, :)));
        end
    end
end

%% мутация
function mutated = mutate(population, rate, lb, ub)
    num_individuals = size(population, 1);
    num_vars = size(population, 2);
    mutated = population;

    for i = 1:num_individuals
        for j = 1:num_vars
            if rand < rate
                % Добавляем случайный шум
                mutated(i, j) = mutated(i, j) + (ub(j) - lb(j)) * 0.1 * randn();

                % Ограничиваем значение в пределах границ
                mutated(i, j) = max(lb(j), min(ub(j), mutated(i, j)));
            end
        end
    end
end


%% уравнения динамики + сопряженные переменные
function dydt = dynamics(t, y, sys)
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