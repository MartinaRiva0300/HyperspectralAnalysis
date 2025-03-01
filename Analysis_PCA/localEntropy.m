function E_loc = localEntropy(numofneigh)
% This function calculates the local entropy defined as the entropy of a 3x3
% grid with the only parameter is the number of neighbours of the central
% element. If the number of neighbors is 8 (the whole matrix is occupied!)
% then the local entropy will be 0
% 
% the length of numofneigh tells me how many neighbours there are
% the value of numofneigh(i) tells me how many neighbours those have
    E_loc = zeros(size(numofneigh));
    for i = 1:length(numofneigh)
        if numofneigh(i) < 8
            p = (1+numofneigh(i))/9;            
            E_loc(i) = -(p*log2(p)+(1-p)*log2(1-p));
            %E_loc(i) = -(p*log2(p));
            %E_loc(i) = acos(2*p - 1)/pi; % alternative entropy calculation
        end
    end
    
end