function [E_tot, dE_tot] = entropyRainfall(mask,idx_sorted)

% E_tot explains how the total entropy of a binary mask of a matrix
% changes by switching on one pixel at a time, as if the pixels were raining
% on a Nan matrix and switching the landing site on from Nan to 1. Which elements are
% switched on and in what order is given by idx_sorted.

deltaE = zeros(1, length(idx_sorted));
E_loc = zeros(1, length(idx_sorted));

numofneigh_mat = nan(size(mask));
for i = 1:length(idx_sorted)
    idx = idx_sorted(i); % the index of the pixel I am looking at
    neigh_idx = getNeighbors(idx, size(mask), 1);  % get the indices of the neighbouring pixels in range=1 
    % I count the number of non Nan neighbors and I assign it to the
    % numofneigh_mat matrix in that index
    numOfNeighOfidx = sum(sum(not(isnan(numofneigh_mat(neigh_idx))))); 
    if numOfNeighOfidx >0
        % I update the value of neighbors because a new neighbour has
        % arrived
        numofneigh_mat(neigh_idx) = numofneigh_mat(neigh_idx)+1;
        % i find the non-Nan neighbors among the pixels that have already landed
        % and save their indices in nonNan_neigh
        nonNan_neigh = neigh_idx(not(isnan(numofneigh_mat(neigh_idx))));
        % I update the local entropy contributions 
        deltaE(i) = sum(localEntropy(numofneigh_mat(nonNan_neigh))-localEntropy(numofneigh_mat(nonNan_neigh)-1));
        numofneigh_mat(idx) = numOfNeighOfidx;
        E_loc(i) = localEntropy(numOfNeighOfidx);
    else % numofneigh_mat_tmp = 0 !
        numofneigh_mat(idx) = numOfNeighOfidx; 
        E_loc(i) = localEntropy(numOfNeighOfidx);
    end
end
% the total entropy is the sum of the previous total entropy plus the new
% local entropy plus the deltaE (changes that the new landing has caused)
E_tot = cumsum(E_loc+deltaE);
dE_tot = E_loc+deltaE;
%E_tot_ratio = E_tot./[1:length(E_tot)];
% Draw the figure
% figure
% plot(E_tot./[1:length(E_tot)])
% xlabel("Number of pixels")
% ylabel ("Total Entropy Ratio")

end


