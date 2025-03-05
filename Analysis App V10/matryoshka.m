function M = matryoshka(parity, coeffs)
% This function takes as input parity, a string which may either be "even"
% or "odd", and a coefficients vector "coeffs" and creates "layered" matrices, like the 
% Russian dolls "matryoshka".
% The core takes the value of the first element of the "weights" vector and
% the elements around the core take the values of the subsequent elements
% of that vector. If the parity is "even" the core consists of the four central
% elements of the matrix.
% EXAMPLE 1: M = matryoshka('odd', [1 0.5]);
%       
% M = [0.5000    0.5000    0.5000;
%      0.5000       1      0.5000;
%      0.5000    0.5000    0.5000]
%
% EXAMPLE 2: M = matryoshka('even', [1 0.5]);
%       
% M = [0.5000    0.5000    0.5000    0.5000;
%      0.5000       1         1      0.5000;
%      0.5000       1         1      0.5000;
%      0.5000    0.5000    0.5000    0.5000]

if strcmp(string(parity), "odd")
    dim = length(coeffs)*2-1;
    center = ceil(dim/2);
    
    M = zeros(dim);
    
    for i=0:floor(dim/2)        
        M(center-i,:) = coeffs(i+1);
        M(center+i,:) = coeffs(i+1);
        M(:,center-i) = coeffs(i+1);
        M(:,center+i) = coeffs(i+1);
    end

elseif strcmp(string(parity), "even")
    
    dim = length(coeffs)*2;
    center = dim/2;

    M = zeros(dim);
    for i = 0:dim/2-1
        M(center-i,:) = coeffs(i+1);
        M(center+1+i,:) = coeffs(i+1);
        M(:,center-i) = coeffs(i+1);
        M(:,center+1+i) = coeffs(i+1);
    end
end

end
