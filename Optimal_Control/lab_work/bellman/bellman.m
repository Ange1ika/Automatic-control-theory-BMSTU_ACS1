clear all
clc

% возможные квантованные состояния и управления
range_x = [0.0, 0.5, 1.0, 1.5];
range_u = [-1.0, -0.5, 0.0, 0.5, 1.0];
n = 2; 

function x_k1 = razn(x_k, u_k)
    a = 0;
    b = 1;
    delta_t = 1;
    x_k1 = (1 + a * delta_t) * x_k + b * delta_t * u_k;
end

function J = fun_J(x_k1, u_k)
    lambda = 2;          % вес управления
    x_target = 0.5;      % цель
    delta_t = 1;
    J = (x_k1 - x_target)^2 + lambda * delta_t * u_k^2;
end


% функция линейной интерполяции для нахождения стоимости 
function cost_interp = interpolate_cost(cost_array, range_x, x_next)
% ищем диапазон, в каких индексах находится x
    idx1 = find(range_x <= x_next, 1, 'last');
    idx2 = find(range_x >= x_next, 1, 'first');

    if isempty(idx1) || isempty(idx2)
        cost_interp = inf; % штрафуем за невозможный переход
    elseif range_x(idx1) == x_next
        cost_interp = cost_array(idx1);
    elseif range_x(idx2) == x_next
        cost_interp = cost_array(idx2);
    else
        x1 = range_x(idx1);
        x2 = range_x(idx2);
        c1 = cost_array(idx1);
        c2 = cost_array(idx2);
        cost_interp = c1 + (c2 - c1) * (x_next - x1) / (x2 - x1);
    end
end

function [u_opt, x_traj, u_traj] = dynamic_prog(range_x, range_u, fun_razn, fun_J, n)
    x_0 = 0.0;
    % инициализируем функцию стоимости 
    C = inf(length(range_x), n+1);
    u_opt = zeros(length(range_x), n);

    % Стоимость при нулевом управлении для x(2)
    for i = 1:length(range_x)
        C(i, n+1) = fun_J(range_x(i), 0); 
    end

    % Итерации по шагам
    for k = n:-1:1
        fprintf('Таблица %d\n', k);
        fprintf('%8s%d%s %8s%d%s %8s%d%s %8s%d%d %8s%d%d %10s%d%s%d%s\n', ...
            'x(', k-1,')', ' u(', k-1,')', 'x(', k,')', 'C', k-1, k, 'J*', k-1, k, 'u*(x(',k-1,'),',k-1,')');
        fprintf('--------------------------------------------------------------------------\n');
        for i = 1:length(range_x)
            min_cost = inf;
            best_u = 0;
            for j = 1:length(range_u)
                u_k = range_u(j);
                x_next = fun_razn(range_x(i), u_k);

                if x_next < min(range_x) || x_next > max(range_x)
                    cost = inf;
                else
                    cost = interpolate_cost(C(:, k+1), range_x, x_next) + fun_J(x_next, u_k);
                end

                if cost < min_cost
                    min_cost = cost;
                    best_u = u_k;
                end

                fprintf(' %9.1f %9.1f %9.1f %12.2f %12.2f %12.1f\n', range_x(i), u_k, x_next, cost, min_cost, best_u);
            end
            C(i, k) = min_cost;
            u_opt(i, k) = best_u;
        end
    end

    % Восстановление оптимальной траектории
    x_traj = zeros(1, n+1);
    u_traj = zeros(1, n);
    x_traj(1) = x_0;

    %fprintf('Восстановление траектории:\n');
    for k = 1:n
        idx = find(abs(range_x - x_traj(k)) < 1e-6, 1);
        if ~isempty(idx)
            u_traj(k) = u_opt(idx, k);
            x_traj(k+1) = fun_razn(x_traj(k), u_traj(k));
            fprintf('Шаг %d: x = %.2f, u = %.2f, x_next = %.2f\n', k-1, x_traj(k), u_traj(k), x_traj(k+1));
        else
            x_traj(k+1) = nan;
            u_traj(k) = nan;
        end
    end
end

[u_opt, x_traj, u_traj] = dynamic_prog(range_x, range_u, @razn, @fun_J, n);

% Построение графиков
figure(Color="white");

yyaxis left;
plot(0:n, x_traj, '-o', 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
ylabel('x^* (Состояние)');
hold on;

yyaxis right;
stairs(0:n-1, u_traj, '-o', 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
ylabel('u^* (Управление)');
xlabel('Шаг');

title('Оптимальная траектория и управление');
legend('Траектория x^*', 'Управление u^*', 'Location', 'Best');
grid on;
grid on;
