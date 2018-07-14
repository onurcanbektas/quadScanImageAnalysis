function [beam_x, beam_y, max_row_pos, max_col_pos] = findCoordinateAxes(image, epsilon_x, epsilon_y)

dim_row = size(image,1);
dim_column = size(image,2);

mean_I_matrix = zeros(dim_row, dim_column);
upper_left_xy = [1, 1];
for i = 1:(dim_row - epsilon_y - 1)
	for j = 1:(dim_column - epsilon_x - 1)
		sum_epsilon = 0;
		for k = 0:epsilon_y
			for l = 0:epsilon_x
				sum_epsilon = sum_epsilon + int32( image(i+k,j+l) );
			end
		end
		mean_I_matrix(i + (epsilon_y./2), j + (epsilon_x./2)) = sum_epsilon;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = mean_I_matrix(:);
total_entry = size(a,1);
max_value = max(max(mean_I_matrix));
%%max_entry_positions = zeros(20,2);

%%[Q,R] = quorem(sym(total_entry), sym(dim_row))
for i = 1:total_entry
	if a(i) >= max_value
        [max_col_pos, max_row_pos] = quorem(sym(i), sym(dim_row));
	end
end

hline(max_row_pos,'r','Analysed Row')
vline(max_col_pos,'b','Analysed Column')

beam_x = image(max_row_pos,:);
beam_y = image(:,max_col_pos)';

end