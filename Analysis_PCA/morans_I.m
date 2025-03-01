function [I_global, I_local, spatial_lag, A_clean] = morans_I(A, weights)

% - I_global is the global Moran's I;
% - I_local is a matrix the same size of A in which in each element the local
% I is stored;
% - B is a matrix the same size of A in which in each element only the
% contribution of the neighbors of that element is stored;
% - A_clean is the same matrix as A-mean(A), except in which each element
% that is completely surrounded by Nan is set to Nan as well, since it
% probably consists of noise.

avgA = mean(A(not(isnan(A)))); % souldn't this be a mean(mean)?
A = A-avgA;

% set the center element of the weights matrix to 0
weights(ceil(size(weights,1)/2), ceil(size(weights,2)/2)) = 0;

% create a Nan frame for the A matrix:
A_framed = nan(size(A,1)+size(weights,1)-1, size(A,2)+size(weights,2)-1);
% insert the A matrix into the frame
A_framed(ceil(size(weights,1)/2) : size(A_framed,1)-floor(size(weights,1)/2),...
   ceil(size(weights,2)/2) : size(A_framed,2)-floor(size(weights,2)/2) ) = A;
clear A

spatial_lag = nan(size(A_framed));

for i = ceil(size(weights,1)/2) : size(A_framed,1)-floor(size(weights,1)/2)
    for j = ceil(size(weights,2)/2) : size(A_framed,2)-floor(size(weights,2)/2)
         if isnan(A_framed(i,j))
            continue
         else
             % create a supporting matrix of the neighbors of the element
             % we are focusing on, which is in the middle
             neighborhood = A_framed( i-floor(size(weights,1)/2) : i+floor(size(weights,1)/2), ...
                   j-floor(size(weights,2)/2) : j+floor(size(weights,2)/2) ); 
             % if all the neighbors are Nan, also the pixel becomes Nan
             if sum(sum(isnan(neighborhood))) == numel(weights)-1                 
                 A_framed(i,j) = nan;
                 continue
             end
            % normalize the weights based on the weights of the non-Nan neighboring elements
            weights_row_normalized = weights./sum(sum(weights(not(isnan(neighborhood)))));
            % bring the nan in the neigh matrix to 0
            neighborhood(isnan(neighborhood)) = 0; 
        end
        spatial_lag(i,j) = sum(sum(neighborhood.*weights_row_normalized));
    end
end

% remove the frames, but maintain the name for memory purposes
A_framed = A_framed(ceil(size(weights,1)/2) : size(A_framed,1)-floor(size(weights,1)/2),...
   ceil(size(weights,2)/2) : size(A_framed,2)-floor(size(weights,2)/2) ) ;
spatial_lag = spatial_lag(ceil(size(weights,1)/2) : size(spatial_lag,1)-floor(size(weights,1)/2),...
   ceil(size(weights,2)/2) : size(spatial_lag,2)-floor(size(weights,2)/2) ) ;

A_clean = A_framed;

%Check whether it's true
varA_clean = var(A_clean(not(isnan(A_clean))));
A_clean = A_clean/sqrt(varA_clean);
spatial_lag = spatial_lag/sqrt(varA_clean);

%% calculation of the local and global Moran's I 
% The formula for I_global states:
%       I_global = sum_i(z_i * sum_j(weights.*neighbors)) / sum_i(z_i.^2)
%
% The formula for I_local states: 
%       I_local = z_i * sum_j(weights.*neighbors) / m2
% where m2 = sum_i(z_i.^2)/N.
% where N is the total number of elements.
%
% It results that I_global = sum_i(I_local)/N

numEl = (sum(sum(not(isnan(A_clean)))));
m2 = (sum( sum((A_clean(not(isnan(A_clean)))).^2)) )/numEl; % isn't this the variance?

% I calculate I_local in two steps to avoid numeric errors amplification
% I start by considering I_local = z_i * (weights.*neighbors) and then I
% will divide by m2
I_local = A_clean.*spatial_lag;

I_global = ( sum(sum(I_local(not(isnan(I_local)))))) / (m2*numEl) ;

I_local = I_local/m2;

%%%%%%%%
% m2 is just the variance of the dataset, which is 1, since I normalize
% it which means that T_global is just the average of the I_locals and
% I_local is the elementwise multiplication of A and the spatial lag


end
