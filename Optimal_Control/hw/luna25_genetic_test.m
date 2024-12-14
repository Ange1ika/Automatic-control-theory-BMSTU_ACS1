function luna25_genetic_general
    % Параметры системы
    sys.tf = 2; 
    sys.epsilon = 1e-4; 
    sys.max_gen = 50; 
    sys.pop_size = 100; 
    sys.mutation_rate = 0.8; 
    sys.crossover_rate = 0.4; 
    sys.n = 2; 

 
    sys.x0 = [1; 2]; 
    sys.xf = [3; 0];

    
    lambda_bounds = [-1,1];

    
    population = initialize_population(sys.pop_size, sys.n, lambda_bounds);

    % ГА
    for gen = 1:sys.max_gen

        % Оценка целевой функции
        fitness = evaluate_population_general(population, sys);

        % Отбор
        selected = tournament_selection(population, fitness, sys.pop_size);

        % Кроссовер
        offspring = crossover(selected, sys.crossover_rate, lambda_bounds);

        offspring = mutate(offspring, sys.mutation_rate, lambda_bounds);
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

    
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
    y0 = [sys.x0; best_lambda(:)];
    [t, sol] = ode45(@(t, y) dynamics_general(t, y), [0 sys.tf], y0, options);
    disp('Полученный xf:');
    disp(sol(:, end));
    visualize_results_general(t, sol);
end

function [t, xf] = simulate_dynamics(sys, lambda)
    y0 = [sys.x0; lambda(:)];
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
    [t, sol] = ode45(@(t, y) dynamics_general(t, y), [0 sys.tf], y0, options);


    xf = sol(:, 1:sys.n)'; 
end

function visualize_results_general(t, sol)
    n = 2;
    figure;
    tiledlayout(n, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:n
        nexttile;
        disp(length(t));
        disp(length(sol(:, i)));
        plot(t, sol(:, i), 'LineWidth', 2);
        grid on;
        title(sprintf('Значение x%d', i), 'FontSize', 12);
        xlabel('Время', 'FontSize', 10);
        ylabel(sprintf('x%d', i), 'FontSize', 10);
    end
end

function population = initialize_population(pop_size, num_params, bounds)
    min_bounds = bounds(1);
    max_bounds = bounds(2);
    population = min_bounds + (max_bounds - min_bounds) * rand(pop_size, num_params);
end

function fitness = evaluate_population_general(population, sys)
    % Оценка невязки
    num_individuals = size(population, 1);
    fitness = zeros(num_individuals, 1);
    for i = 1:num_individuals
        lambda = population(i, :)';
        [~, xf] = simulate_dynamics(sys, lambda);
        fitness(i) = norm(xf - sys.xf);
    end
end

function selected = tournament_selection(population, fitness, pop_size)
    num_individuals = size(population, 1);
    selected = zeros(pop_size, size(population, 2));
    for i = 1:pop_size
        idx = randi(num_individuals, 2, 1);
        if fitness(idx(1)) < fitness(idx(2))
            selected(i, :) = population(idx(1), :);
        else
            selected(i, :) = population(idx(2), :);
        end
    end
end

function offspring = crossover(population, rate, bounds)
    num_individuals = size(population, 1);
    num_vars = size(population, 2);
    offspring = population;
    for i = 1:2:num_individuals-1
        if rand < rate
            mask = rand(1, num_vars) > 0.5;
            parent1 = population(i, :);
            parent2 = population(i+1, :);
            offspring(i, :) = parent1 .* mask + parent2 .* ~mask;
            offspring(i+1, :) = parent2 .* mask + parent1 .* ~mask;
        end
    end
end

function mutated = mutate(population, rate, bounds)
    num_individuals = size(population, 1);
    num_vars = size(population, 2);
    mutated = population;
    for i = 1:num_individuals
        for j = 1:num_vars
            if rand < rate
                mutated(i, j) = mutated(i, j) + (bounds(2) - bounds(1)) * 0.2 * randn();
                mutated(i, j) = max(bounds(1), min(bounds(2), mutated(i, j)));
            end
        end
    end
end

function dydt = dynamics_general(~, y)
    % Уравнения сопряжённых переменных
    % x1 = y(1);         
    x2 = y(2);         
    lambda1 = y(3); 
    lambda2 = y(4); 
   
    dx1 = x2;                
    dx2 = -lambda2;          
    dlambda1 = 0;            
    dlambda2 = -lambda1;     
    
    dydt = [dx1; dx2; dlambda1; dlambda2];
end


