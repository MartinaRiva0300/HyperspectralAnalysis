function [r,C]=rcalc(kinit1,c0,tf,cexp,choose,O,R,Ccol)

% fit differential equations
C=cell(size(c0,1),1);
Residual=[];
for i=1:size(c0,1)
    [t_ode,C{i,1}]=ode45('equations',tf{i,1},c0(i,:),odeset('RelTol',1e-2,'AbsTol',1e-3),O,R,abs(kinit1),choose);
    cinput=cexp{i,1};
    Cfit=C{i,1};
    Rin=cinput(:,logical(Ccol))-Cfit(:,logical(Ccol));
    Residual=[Residual;Rin];
end
% Global residual to be minimised
r=Residual(:);
