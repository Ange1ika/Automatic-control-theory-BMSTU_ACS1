function lunar_landing_shooting

    % Задание начальных условий
    t0 = 0;
    tf = 500;
    x0 = [0; 18000; 0; 1690; 1200]; % Начальные состояния
    xf = [NaN; 50; 2; 5; NaN];       % Терминальные условия
    mu = 4902.8;                    
    P = 450;                        
    K = 1000;
    eps = 0.1;                      % Допуск по погрешности пристрелки
    max_iter = 10000;                  

    lambda0 = [0; 0; 0; 0; 0];  % Начальное предположение для λ(0)
    alpha_prev = 0;              % Начальное значение угла управления
    prev_error = zeros(size(lambda0)); % Инициализация предыдущей ошибки для корректировки lambda

    for iter = 1:max_iter
      
        % Задание начальных условий для интеграции
        init_conditions = [x0; lambda0; alpha_prev];
        
        % Интеграция динамики
        [t, x, lambda, alpha] = integrate_dynamics(t0, tf, init_conditions, mu, P, K, alpha_prev);
        
        % Проверка терминальных условий
        terminal_errors = compute_terminal_errors(x(:, end), xf);
        if norm(terminal_errors) < eps
            fprintf('Терминальные условия выполнены за %d итераций\n', iter);
            break;
        end
        
        % Обновление λ(0) с использованием ошибки, игнорируя NaN в терминальных условиях
        lambda0 = lambda0 + eps * terminal_errors - 0.5 * (terminal_errors - prev_error);
        %disp(lambda0);
        prev_error = terminal_errors; % Обновление предыдущей ошибки
        %disp(alpha);
        alpha_prev = mean(alpha);

        fprintf('Итерация %d: Ошибка терминальных условий = %f\n', iter, norm(terminal_errors));
    end
    
    if iter == max_iter
        warning('Достигнуто максимальное число итераций. Результат может быть неточным.');
    end
    
    plot_results(t, x, lambda, alpha);
end

%% Решение ДУ

function [t, x, lambda, alpha] = integrate_dynamics(t0, tf, init_conditions, mu, P, K, alpha_prev)
    % Задаем параметры интегратора
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    
    [t, X] = ode45(@(t, X) dynamics(t, X, mu, P, K, alpha_prev), [t0 tf], init_conditions, options);
    
    x = X(:, 1:5)';
    %disp(x);
    lambda = X(:, 6:10)';
    alpha = X(:, 11)';
end

%% Расчёт производных Гамильтониана по переменным состояния для сопряжённых уравнений

function dlambdadt = jacobian_H(x, lambda, alpha, mu, P, K)
    
    phi = atan2(x(1), x(2));
    r = sqrt(x(1)^2 + x(2)^2);

    dlambda1 = 2 * x(1) + lambda(3) * (2 * x(1) * mu * sin(phi) / r^2) + lambda(4) * (2 * x(1) * mu * cos(phi) / r^2);
    dlambda2 = 2 * x(2) + lambda(3) * (2 * x(2) * mu * sin(phi) / r^2) + lambda(4) * (2 * x(2) * mu * cos(phi) / r^2);
    dlambda3 = 2 * x(3) + lambda(1);
    dlambda4 = 2 * x(4) + lambda(2);
    dlambda5 = 2 * x(5) + lambda(3) * P * (cos(alpha) * cos(phi) - sin(alpha) * sin(phi)) / x(5)^2 - P * lambda(4) * (sin(alpha) * cos(phi) + cos(alpha) * sin(phi)) / x(5)^2;

    dlambdadt = [dlambda1; dlambda2; dlambda3; dlambda4; dlambda5];
end

%% Расчёт уравнений динамики

function dXdt = dynamics(t, X, mu, P, K, alpha_prev)
    x = X(1:5);
    lambda = X(6:10);
    
    phi = atan2(x(1), x(2));
    alpha = atan2(lambda(4) * cos(phi) - lambda(3) * sin(phi), cos(phi) + lambda(4) * sin(phi) / lambda(3));
    alpha = (alpha + alpha_prev) / 2;
    
    dxdt = [x(3);
            x(4);
            -mu/(x(1)^2 + x(2)^2) * sin(phi) - (P * cos(alpha) * cos(phi)) / x(5) + (P * sin(alpha) * sin(phi)) / x(5);
            -mu/(x(1)^2 + x(2)^2) * cos(phi) + (P * sin(alpha) * cos(phi)) / x(5) + (P * cos(alpha) * sin(phi)) / x(5);
            -K * P];
    
    dlambdadt = -jacobian_H(x, lambda, alpha, mu, P, K);
    
    dXdt = [dxdt; dlambdadt; alpha]; 
end

%% Расчёт ошибки по терминальным условиям

function terminal_errors = compute_terminal_errors(xf_actual, xf_desired)
    terminal_errors = xf_actual - xf_desired;
    terminal_errors(isnan(xf_desired)) = 0; % Не учитываем NaN условия
end

%% Отрисовка графиков

function plot_results(t, x, lambda, alpha)
    figure;
    subplot(3, 1, 1);
    plot(t, x(1, :), t, x(2, :), t, x(3, :), t, x(4, :), t, x(5, :));
    legend('x', 'y', 'V_x', 'V_y', 'm');
    title('Траектория движения и состояния');
    % Расчёт производных Гамильтониана по переменным состояния для сопряжённых уравнений
    subplot(3, 1, 2);
    plot(t, lambda(1, :), t, lambda(2, :), t, lambda(3, :), t, lambda(4, :), t, lambda(5, :));
    legend('\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4', '\lambda_5');
    title('Сопряженные переменные');
    
    subplot(3, 1, 3);
    plot(t, alpha);
    legend('\alpha');
    title('Угол управления');
end
