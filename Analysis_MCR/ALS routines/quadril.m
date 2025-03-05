function [c,u1,u2,u3]=quadril(c,ne1,ne2,ne3,ishape)
% function [c,u1,u2,u3]=quadril(c,ne1,ne2,ne3,ishape)
% forces the shapes of the profiles (conc) to follow the quadrilinear model
% It estimates the profiles in the three modes concatenated in the bilinear
% profile mode: c(nst,1) --> u1(ne1,1),u2(ne2,1),u3(ne3,1)
% nst=ne1*ne2*ne3
% c input-output unfolded bilinear profile

[nst,dum]=size(c);

% ********************************************************
% change the order of the dimensions to do the appropriate
% refolding and svd; ne1 is the dimensions of mode 1, which 
% is the row size of the individual augmented matrices 
%
dummy=ne1;ne1=ne3;ne3=dummy;
%
% ********************************************************


% first check mode dimensions are equal to length of profile

if ne1*ne2*ne3 ~= nst,
    error('mode dimensions do not fit the lenght of the profile'),
end

% inizialization
ns1=nst/ne1;
ni=1;
nf=ns1;

% folds the unfolded profiles in a matrix
% sinchronization ishape=1, no synchronization ishape=2

for j=1:ne1,
    cf1(1:ns1,j)=c(ni:nf,1);
    ni=nf+1;
    nf=ni+ns1-1;
end

if ishape==2,
    y=peakshift(cf1);
else
    y=cf1;
end

% **********************************************************
% svd of the first folded profiles matrix
% **********************************************************

[u,s,v]=svd(y,0);

% estimates the profiles in the first mode
u1=s(1,1).*v(:,1)';
cf1c=u(:,1)*u1;
rs=cf1-cf1c;
un=sum(sum(rs*rs'));
sigma=sqrt(un/nst);
disp('std. dev. res. between the estimations');disp(sigma);
disp('singular values ');disp([s(1,1);s(2,2)]);
% pause

% ************************************************************
% calculates the second folded profiles matrix
% **************************************************************

% [nst,dum]=size(u(:,1));
ns2=ns1/ne2;
ni=1;
nf=ns2;

% folds the unfolded profiles in a matrix
% assumes there is sinchronization

for j=1:ne2,
    cf2(1:ns2,j)=u(ni:nf,1);
    ni=nf+1;
    nf=ni+ns2-1;
end
y=cf2;

% **********************************************************
% new svd of the profile matrix
% **********************************************************

[u,s,v]=svd(y,0);

% estimates the total concentrations
u2=s(1,1).*v(:,1)';
u3=u(:,1);
cf2c=u3(:,1)*u2;
rs=cf2-cf2c;
un=sum(sum(rs*rs'));
sigma=sqrt(un/ns1);
disp('std. dev. res. between the estimations');disp(sigma);
disp('singular values ');disp([s(1,1);s(2,2)]);
% pause

% *************************************************
% recalculates the new long refolded profile
% *************************************************

y=cf2c;
ni=1;
nf=ns2;
for j=1:ne2,
    cf2r(ni:nf,1)=y(1:ns2,j);
    ni=nf+1;
    nf=ni+ns2-1;
end

y=cf2r*u1;

% refold the profiles
ni=1;
nf=ns1;
for j=1:ne1,
    cf1r(ni:nf,1)=y(1:ns1,j);
    ni=nf+1;
    nf=ni+ns1-1;
end


% calculates the differences
rs=cf1r-c;
u=sum(sum(rs*rs'));
sigma=sqrt(u/nst);
disp('std. dev. res. between the estimations');disp(sigma);
%disp('singular values ');disp([s(1,1);s(2,2)]);

%subplot(2,1,1),plot(cf1r)
%subplot(2,1,2),plot(c,'r')
%pause(1)
% final assigment

c=cf1r;

end

