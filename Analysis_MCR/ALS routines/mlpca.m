function [U,S,V,SOBJ,ErrFlag] = mlsmall(X,stdX,p,convlim,maxiter);
%
%                 MLPCA.M   v. 4.0
%
% This function performs maximum likelihood principal components
% analysis assuming uncorrelated measurement errors.  The variables
% passed to the function are:
%
%   X    is the mxn matrix of observations (measurements).
%
%   stdX is the mxn matrix of standard deviations associated with
%        the observations in xobs.  For missing measurements, stdX
%        should be set to zero.
%
%   p    is the dimensionality of the model (p<=n and p<=m).
%
% The parameters returned are:
%
%   U,S,V  are pseudo-svd parameters (mxp, pxp, and nxp).  The
%          maximum likelihood estimates are given by:
%                   XML=U*S*V'
%          NOTE: When using the results of MLPCA to project new
%          measurements onto the model (e.g. in PCR), orthogonal
%          projection is not recommended, since this defeats the
%          purpose of maximum likelihood estimation.  Instead the
%          projection should be done by:
%               xml=U*inv(U'*Q*U)*U'*Q*x or
%               xml=V*inv(V'*Q*V)*V'*Q*x
%          depending on whether the measurements are in the row
%          or column space of the original vector.  Q is the
%          inverse of the mxm or nxn diagonal matrix of measurement
%          variances.
%
%   SOBJ is the value of the objective function for the best model.
%        For correct error estimates, this should follow a
%        chi-squared distribution with (m-p)*(n-p) degrees of
%        freedom.
%
%   ErrFlag indicates the termination conditions of the function;
%             0 = normal termination (convergence)
%             1 = maximum number of iterations exceeded
%
% The function can also produce a file - mlpca.mat - if appropriate
% lines are activated as indicated in the code.  This can be used to
% follow convergence if desired.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization
%
% convlim=1e-10;             % convergence limit
% maxiter=200000;            % maximum no. of iterations
XX=X;                      % XX is used for calculations
varX=(stdX.^2);            % convert s.d.'s to variances
[i,j] = find(varX==0);     % find zero errors and convert to large
errmax = max(max(varX));   % errors for missing data
for k=1:length(i);
   varX(i(k),j(k)) = 1e+10*errmax;
end
n=length(XX(1,:));         % the number of columns
%
% Generate initial estimates assuming homoscedastic errors.
%
[U,S,V]=svd(XX,0);         % Decompose adjusted matrix
U0=U(:,1:p);               % Truncate solution to rank p
count=0;                   % Loop counter
Sold=0;                    % Holds last value of objective function
ErrFlag=-1;                % Loop flag
%
% Loop for alternating regression
%
while ErrFlag<0;
   count=count+1;          % Increment loop counter
%
% Evaluate objective function
%
   Sobj=0;                             % Initialize sum      
   MLX=zeros(size(XX));                % and maximum likelihood estimates
   for i=1:n                           % Loop for each column of XX
      Q=sparse(diag(varX(:,i).^(-1))); % Inverse of error covariance matrix
      F=inv(U0'*Q*U0);                    % Intermediate calculation
      MLX(:,i)=U0*(F*(U0'*(Q*XX(:,i))));  % Max. likelihood estimates
      dx=XX(:,i)-MLX(:,i);                % Residual vector
      Sobj=Sobj+dx'*Q*dx;                 % Update objective function
   end
%
% This section for diagnostics only and can be commented out.  "Ssave"
% can be plotted to follow convergence.
%
%   Ssave(count)=Sobj;
%   save mlpca;
%
% End diagnostics
%
% Check for convergence or excessive iterations
%
   if rem(count,2)==1                   % Check on odd iterations only
      if (abs(Sold-Sobj)/Sobj)<convlim  % Convergence criterion
         ErrFlag=0;
      elseif count>maxiter              % Excessive iterations?
         ErrFlag=1;
      end
   end
%
% Now flip matrices for alternating regression
%
   if ErrFlag<0                         % Only do this part if not done
      Sold=Sobj;                        % Save most recent obj. function
      [U,S,V]=svd(MLX,0);               % Decompose ML values
      XX=XX';                           % Flip matrix
      varX=varX';                       % and the variances
      n=length(XX(1,:));                % Adjust no. of columns
      U0=V(:,1:p);                      % V becomes U in for transpose
   end
end
%
% All done.  Clean up and go home.
%
[U,S,V]=svd(MLX,0);
U=U(:,1:p);
S=S(1:p,1:p);
V=V(:,1:p);
SOBJ=Sobj;



 
 
 
 
 

