function [cnew,tnew]=tuck2(c,nmat,ishape,ituck)
% function [cnew,t]=tuck(c,nmat,ishape)
% c input-output unfolded profiles
% nmat number of simultaneous matrices/experiments considered
% ishape controls THREE-WAY TYPE OF MODEL CONSTRAINTS 
% Non trilinear, only bilinear (0)
% Trilinear, equal shape and synchronization (all species (1)
% Trilinear without synchronization (all species), (2)
% Trilinear and synchronization for only some species, (3)
% Dif. nr. of components in the three modes, TUCKER3 models, (4)

% nrt nr. of total rows
% ncommon nr. of profiles interacting 
[nrt,ncomp]=size(c);

% nr. of rows for individual matrix
nr=nrt/nmat;
% ishape = 4
% refolds the unfolded vector profiles in a single matrix

if ishape==4
    cp=[]; 
    cp1=zeros(nr,1);

    if ituck==1,
        for j1=1:ncomp
            ni=1;
            nf=nr;
            for j2=1:nmat,
                cp1=c(ni:nf,j1);
                cp=[cp,cp1]
                ni=nf+1;
                nf=ni+nr-1;
            end
        end
    end
    if ituck==3,
        for j1=1:ncomp
            ni=1;
            nf=nr;
            cp2=[];
            for j2=1:nmat,
                cp1=c(ni:nf,j1);
                cp2=[cp2,cp1];
                ni=nf+1;
                nf=ni+nr-1;
            end
            cp=[cp;cp2];
        end
    end
y=cp;

[u,s,v]=svd(y,'econ');

t=s(1,1).*v(:,1);

k=0;
if ituck==1,
    for i=1:ncomp,
        for j=1:nmat,
            k=k+1;tnew(j,i)=t(k);
        end,
    end
end

if ituck==3,tnew=t;end
    
%  size(t)
%  display(tnew)
%  size(tnew)
%  size(u)
%  size(v)
%  ituck
 
% calculates the new profiles using kronecker product

if ituck==1,cnew=kron(tnew,u(:,1));end

if ituck==3,
    unew=[];ni=1;nf=nr;
    for j1=1:ncomp
        unew(:,j1)=u(ni:nf,1);
        ni=nf+1;nf=ni+nr-1;
    end
    cnew=kron(tnew,unew);
end

% size(cnew)

l=svd(cnew);
if l(ncomp)<10-5, disp('rank deficient model'),pause,end
figure; plot(cnew)
pause(1)
close
end
