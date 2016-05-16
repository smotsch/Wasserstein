function [A, B] = AdjustSize(A,B)
% This function takes in two vectors and returns two vectors of the same
% length by adding NaN to the shorter vector.

% Author: Shane Lubold (shane.lubold@asu.edu)

    if length(A) == max(length(A),length(B))
        B = [B; NaN(1,abs(length(A) - length(B)))'];
    else
        A = [A; NaN(1,abs(length(A) - length(B)))'];
    end
end

