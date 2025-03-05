function [r2t2,loft]=loftril(dt,u,v,t)
% function [r2t2,loft]=loftril(dt,u,v,t)
% dt should be entered as a cube
% u(nr,n), v(nc, n),t(nmat, n) loading matrices in the three modes 
% entered always as column matrices
% nr n.ro3ws, nc nr. columns, nmat nr matrices or slices, n nr. components

[nr,nc,nt]=size(dt);
[nc,n]=size(v);

if nt==1, 
    for i=1:nr
        for j=1:nc
            dtc(i,j)=0;
            dtc(i,j)=dtc(i,j)+u(i,:)*v(j,:)';
            res(i,j)=dt(i,j)-dtc(i,j);
        end
    end
    
% calculate the residuals and total sum of suares
sumdt=sum(sum(dt.*dt));
sumdtc=sum(sum(dtc.*dtc));
sumres=sum(sum(res.*res));
else
    
% calculated repreduced data cube and residuals
    for i=1:nr
        for j=1:nc
            for k=1:nt
                dtc(i,j,k)=0;
                for l=1:n
                    dtc(i,j,k)=dtc(i,j,k)+u(i,l)*v(j,l)*t(k,l);
                end
                res(i,j,k)=dt(i,j,k)-dtc(i,j,k);
            end
        end
    end

% calculate the residuals and total sum of suares
sumdt=sum(sum(sum(dt.*dt)));
sumdtc=sum(sum(sum(dtc.*dtc)));
sumres=sum(sum(sum(res.*res)));

%pause

end

disp(['sstot,sscalc and ssres = ',num2str([sumdt,sumdtc,sumres])]);
loft=(sqrt(sumres/sumdt))*100;
r2t1=((sumdtc)/sumdt)*100;
r2t2=((sumdt-sumres)/sumdt)*100;
disp(['lof (%) = ',num2str(loft)]);
disp(['fit (%) = ',num2str(100-loft)]);
% disp(['R square is = ',num2str(r2)]);
% fprintf(1,'Rsquare is %15.12f\n',r2t1',r2t2)
disp(['expl.var (unique)% = ',num2str(r2t1)]);
disp(['R2 % = ',num2str(r2t2)]);

disp([ ])

end
