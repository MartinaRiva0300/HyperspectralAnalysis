function varargout = nn_shuffle(M)
% This function shuffles non-Nan elements in a matrix and leaves the
% positions of the Nan elements unchanged.
% 1st OUTPUT: it's the shuffled matrix
% 2nd OUTPUT: 2 column vector in which the first column is the non-Nan
% indices of the original matrix and the second column is those indices
% randomized. 
% EXAMPLE: If the first element of the first column is 2 and that of
% the second column is 5, then that means that the first non-Nan element
% of the original matrix is in position 2. In that same position, in the
% shuffled matrix we find the element that was in position 5 in the
% original.

idx_matrix = reshape(1:size(M,1)*size(M,2), size(M));

ord_idx = idx_matrix(not(isnan(M))); 
tmp =randperm(length(ord_idx));
rnd_idx = ord_idx(tmp);

M_shuffled = nan(size(M));
M_shuffled(ord_idx) = M(rnd_idx);

for k = 1:nargout
    if k == 1
    varargout{k} = M_shuffled;
    end
    if k == 2
    varargout{k} = cat(2, ord_idx, rnd_idx);
    end
end