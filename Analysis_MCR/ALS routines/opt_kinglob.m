function [kout,conckout,ssq_nglm,J_nglm]=opt_kinglob(cexp,R,O,c0,kinit,finaltime,choose,Ccol)

% [k,c]=opt_kin(cexp,m,n,vec,c0,kinit,tf,choose,nmod)
%
% Determines the optimal rate constants k and concentration profiles
% 
% cexp   - experimental concentration profiles
% R      - matrix of stoichiometric coefficients
% O      - matrix of order of reactions
% c0     - vector with initial concentration of species
% kinit  - vector with initial estimate of rate constants
% finaltime - time axis
% choose - vector with species to be considered for each reaction
%
% k      - optimal values of rate constants
% c      - optimal concentration profiles (determined through k)

%*********************************
disp('opt_kin')
%*********************************

% Matrix kinit - Initial estimates of kinetic constants ****
kinit1=kinit';
kinit1=kinit1(find(kinit1));

% Optimization of k - Rate constants ***********************

[kout,conckout,ssq_nglm,J_nglm]=nglmglob(@rcalcglob,kinit1,c0,finaltime,cexp,choose,O,R,Ccol);
kout=abs(kout);

% Determination of c - Concentration profiles **************
C=cell(size(c0,1),1);
for i=1:size(c0,1)
    [t,C{i,1}]=ode45('equations',finaltime{i,1},c0(i,:),odeset('RelTol',1e-2,'AbsTol',1e-3),O,R,kout,choose);
end

