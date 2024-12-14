possible_moves = zeros(3, 3);
%possible_moves(:, 1) = ones(3, 1);

%disp(possible_moves)

%disp(possible_moves(:, 3))
num_rows_map = 3;
num_cols_map = 3;

padded_map = ones(num_rows_map + 2, num_cols_map + 2);
padded_map(2:end-1, 2:end-1) = possible_moves;

disp(padded_map)