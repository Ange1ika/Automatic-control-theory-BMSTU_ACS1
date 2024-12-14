%% для 2

syms p11 p12 p22

P = [p11 p12; p12 p22];

A = [0 1; 2 -1];
B = [0; 1];
Q = [2 0; 0 1];     % Весовая матрица ошибок регулирования
R = [1/2];          % Весовая матрица управления

x0 = [-4 4];
% lqr - выдается матрица оптимальных коэф усиления
% P - решение Риккарти 
% - P*A - A'*P - Q + P*B*inv(R)*B'*P


[K, P, EV] = lqr(A, B, Q, R)

C = [1 1];
D = [1];
Bin = [0; 0]; % фиктивная матрица для initial

h = 0.01;
t_min = 0;
t_max = 10;
% t = (t_min : h : t_max);

% Реакция системы на начальное воздействие в пространстве состояний
[Y, X, t] = initial(A-B*K, Bin, C, D, x0, t_max);

u = -K * X';
u_1 = -inv(R) * B' * P * X';


 


% labels = {'x (Положение по оси X)', 'y (Положение по оси Y)'};
figure(color = 'white')
nexttile;
plot(t, X(:, 1), 'LineWidth', 2, 'Color', 'r'); % Разные цвета для каждой линии
grid on;
hold on;
plot(t, X(:, 2), 'LineWidth', 2, 'Color', ' b');

title('x (Положение по оси X)', 'FontSize', 14, 'FontWeight', 'bold'); % Увеличиваем шрифт заголовка
xlabel('Время, с', 'FontSize', 12); % Подписываем оси с единицами измерения
%ylabel(labels{i}, 'FontSize', 12);  % Подписываем оси с единицами измерения
ax = gca;
ax.FontSize = 12; % Увеличиваем шрифт на осях
ax.LineWidth = 1.5; % Утолщаем линии осей
ax.GridLineStyle = '--'; % Делаем сетку пунктирной



% labels = {'x (Положение по оси X)', 'y (Положение по оси Y)'};
figure(color = 'white')
nexttile;
plot(t, u_1, 'LineWidth', 2, 'Color', 'r'); % Разные цвета для каждой линии
grid on;
%hold on;
%plot(t, X(:, 2), 'LineWidth', 2, 'Color', 'p');

title('Управление', 'FontSize', 14, 'FontWeight', 'bold'); % Увеличиваем шрифт заголовка
xlabel('Время, с', 'FontSize', 12); % Подписываем оси с единицами измерения
%ylabel(labels{i}, 'FontSize', 12);  % Подписываем оси с единицами измерения
ax = gca;
ax.FontSize = 12; % Увеличиваем шрифт на осях
ax.LineWidth = 1.5; % Утолщаем линии осей
ax.GridLineStyle = '--'; % Делаем сетку пунктирной