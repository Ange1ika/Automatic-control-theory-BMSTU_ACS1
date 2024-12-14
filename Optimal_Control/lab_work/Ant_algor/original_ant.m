clc
clear all

% 1.1 Создаем модель
stAntOpts = struct();

% Задание параметров модели
stAntOpts.NumAnts = 10;
stAntOpts.MaxIter = 30;
stAntOpts.Alpha = 2; 
% альфа - коэффициент феромона, при 0 будем ориентироваться только на
% кратчайший путь 

stAntOpts.Beta = 1;
% бета - коэффициент расстояния, при 0 будем
% ориентироваться только на оставляемый запах

stAntOpts.PheromoneEvaporation = 0.5;


%% 1.2 создание бинарной матрицы с препятствиями

function grid = create_map(size_x, size_y, barrier_coord)
    grid = zeros(size_x, size_y);

    for i = 1:(length(barrier_coord))
        x = barrier_coord(i, 1);
        y = barrier_coord(i, 2);
        if x >= 1 && x <= size_x && y >= 1 && y <= size_y
            grid(x, y) = 1;
        else
            warning("Координаты (%d, %d) выходят за границы")
        end
    end
end


%% вычисление пути от заданной начальной точки к 
% заданной конечной точке по карте с 
% использованием муравьиного алгоритма 
% 1.3
function path = ant_colony_optimization(start_point, end_point, ...
    map, stAntOpts)

    [size_x, size_y] = size(map);

    % инициализация начальных значений феромонов на карте
    matrix_pheromones = 10 * ones(size_x, size_y);
    path = [];
    
    best_path_len = Inf;

    for i = 1:stAntOpts.MaxIter
        ant_paths = cell(1, stAntOpts.NumAnts);
    
        % вычисляем путь каждого муравья
        for k = 1 : stAntOpts.NumAnts
            curr_point = start_point;
            ant_path = curr_point;

            while ~isequal(curr_point, end_point)
                % Определение след положения
                possible_moves = find_possible_moves(curr_point, ...
                    end_point, map, ant_path);

                if isempty(possible_moves)
                    ant_path = start_point; % Сброс пути
                    break
                end

                probs = calculate_probabilities(possible_moves, ...
                    matrix_pheromones, curr_point, end_point, stAntOpts);
        
                next_move = choose_next_move(probs, possible_moves);
                % обновление текущей точки
                curr_point = next_move;
                ant_path = [ant_path; curr_point];
            end
            % Обновление феромонов после одного муравья
            matrix_pheromones = update_pheromones(matrix_pheromones, ...
                ant_path, end_point, stAntOpts)
            % сохраняем путь муравья в ячейку
            ant_paths{k} = ant_path;

            % Обновление лучшего пути 
            if (size(ant_paths, 1) < best_path_len ...
                    && isequal(ant_path(end, :), end_point))
                path = ant_path;
                best_path_len = size(ant_path, 1);
            end
        end

        % Обновление феромонов после итерации
        matrix_pheromones = update_pheromones_all(matrix_pheromones, ...
            ant_paths, end_point, stAntOpts);

        fprintf("Итерация  %d: длина лучшего пути = %d\n", i, best_path_len); 
    end
end

   



%% 1.4 нахождение возможных перемещений
function possible_moves = find_possible_moves(current_point, ...
    end_point, map, path)
    
    [num_rows_map, num_cols_map] = size(map);
    
  
    x_cur = current_point(1);
    y_cur = current_point(2);
    
    % Паддинг карты с единицами
    padded_map = ones(num_rows_map + 2, num_cols_map + 2);
    padded_map(2:end-1, 2:end-1) = map;

    
    x_cur = x_cur + 1;
    y_cur = y_cur + 1;
    
   
    possible_moves = [];
    
    % смещения для 8 направлений
    directions = [-1, -1; -1, 0; -1, 1; 0, -1; 0, 1; 1, -1; 1, 0; 1, 1];
    
 
    for i = 1:size(directions, 1)
        new_x = x_cur + directions(i, 1);
        new_y = y_cur + directions(i, 2);
        
        
        if padded_map(new_x, new_y) == 0
            orig_x = new_x - 1;
            orig_y = new_y - 1;
            
            % Проверяем, что новое положение не является частью 
            % уже пройденного пути
            if ~ismember([orig_x, orig_y], path, 'rows')
                possible_moves = [possible_moves; orig_x, orig_y];
            end
        end
    end
end



%% 1.5 расчёт вероятностей для перемещений
function probabilities = calculate_probabilities(possible_moves, pheromones, ...
    current_point, end_point, stAntOpts)

    eps = 1e-3;
    
    x_cur = current_point(1);
    y_cur = current_point(2);

    num_moves = size(possible_moves, 1);
    %disp(num_moves);
    probabilities = zeros(num_moves, 1);

    for i = 1:num_moves
        move = possible_moves(i, :);
        x_move = move(1);
        y_move = move(2);

        tau = pheromones(x_move, y_move);
        dist = sqrt((end_point(1)- x_move)^2 + (end_point(2)-y_move)^2);
        eta = 1/(dist + eps);

        probabilities(i, :) = tau^stAntOpts.Alpha * eta^stAntOpts.Beta;
    end

    probabilities = probabilities/sum(probabilities);

end

%% 1.6 расчёт следующего перемещения
function next_move = choose_next_move(probabilities, possible_moves)
   
    % Генерируем случайное число от 0 до 1
    r = rand();
    
    % превращаем вероятности в отрезок от 0 до 1
    % (например prob = [0.1, 0.5, 0.4]
    % cum = cumsum(prob) = [0-0.1, 0.1-0.6, 0.6-1])
    cum_probabilities = cumsum(probabilities);
    
    % Определяем индекс следующего перемещения
    index = find(cum_probabilities >= r, 1);
    
    % Координаты следующего перемещения
    next_move = possible_moves(index, :);
    if isempty(index)
        error('Не удалось выбрать следующее перемещение: все вероятности равны нулю.');
    end
end

%% 1.7 обновление феромонов локально для одного внутри цикла
function pheromones = update_pheromones(pheromones, path, end_point, ...
    stAntOpts)
    
    p = stAntOpts.PheromoneEvaporation;
    delta_tau = 1/length(path);

    % Испарение
    pheromones = (1 - p) * pheromones;


    %Распространение — добавление новых феромонов на основе пройденных путей.

    if end_point == path(end, :)

        % Добавление феромонов по пути
        % можно сделать с усилением насыщенности феромона, но тогда
        % вероятен сход к локальному оптимуму

        for i = 1:size(path, 1)
            x = path(i, 1);
            y = path(i, 2);
            pheromones(x, y) = pheromones(x, y) + delta_tau;
        end

    else
        for i = 1:size(path, 1)
            x = path(i, 1);
            y = path(i, 2);
            pheromones(x, y) = pheromones(x, y) - 3;
            pheromones(pheromones < 0) = 0;
        end
    end
end

%% 1.8 обновление феромонов глобально

function pheromones = update_pheromones_all(pheromones, ants_paths, ...
    end_point, stAntOpts)

    p = stAntOpts.PheromoneEvaporation;
    pheromones = (1 - p) * pheromones;

    % Распространение феромонов по путям всех муравьев
    for k = 1:length(ants_paths)
        path = ants_paths{k}; 

       
        if isequal(path(end, :), end_point)
            delta_tau = 1 / size(path, 1); % Усиление феромонов пропорционально длине пути

       
            for i = 1:size(path, 1)
                x = path(i, 1);
                y = path(i, 2);
                pheromones(x, y) = pheromones(x, y) + delta_tau;
            end
        else
        
            for i = 1:size(path, 1)
                x = path(i, 1);
                y = path(i, 2);
                pheromones(x, y) = pheromones(x, y) - 3; 
            end
        end
    end

    pheromones(pheromones < 0) = 0;
end

%% Поиск оптимального пути с использованием муравьиного алгоритма

rng(2696);

size_x = randi([15, 30]);
size_y = randi([15, 30]);

n_obstacles = randi([30, 50]);

barrier_coord_x = randi([1, size_x], n_obstacles, 1);
barrier_coord_y = randi([1, size_y], n_obstacles, 1);
barrier_coord = [barrier_coord_x, barrier_coord_y];

start_point = [1, 1];
end_point = [size_x, size_y];

map_with_obst = create_map(size_x, size_y, barrier_coord);

% ant_colony_optim

path = ant_colony_optimization(start_point, end_point, ...
    map_with_obst, stAntOpts);
plot_map(map_with_obst, path);


% Построение графика бинарной карты с препятствиями и оптимальным путем
    % map  - бинарная карта с препятствиями
    % path - массив с координатами оптимального пути

    function plot_map(map, path)

    [size_y, size_x] = size(map);

    % Добавление на карту с препятствиями оптимального пути
    for i = 1:length(path)
        map(path(i,1), path(i,2)) = 2;   
    end
    
    % Добавление пустой строки и столбца для корректного вывода
    a = zeros(size_y, 1);
    b = zeros(1, size_x + 1);
    
    map = [map a];
    map = [map; b];
    
    
    h = figure(); 
    clf;

        pcolor(map);
        my_map = [1 1 1; 0 0 0; 1 0 0];
        colormap(my_map); % Установка цветовой карты
    
        view(0, -90);
    
        xlabel('x','FontName','TimesNewRoman','FontSize',12);
        ylabel('y','FontName','TimesNewRoman','FontSize',12);
   
end

