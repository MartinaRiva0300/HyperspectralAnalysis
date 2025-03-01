%program that, after the user has loaded the Spectral Hypercube, asks the
%user to input a certain number of frequencies and to choose the ROI from
%the image. This program plots all the complex values in the Gauss plane.
%Each frequency is shown in a different color. The arrows indicate the
%means. The circumferences indicate different type of standard deviation.

MAT_HYP=0; %MAT_HYP=1 load .mat hypercube

if MAT_HYP==1
    dir0=('C:\Users\HARDi\Desktop\Tesi');
    [filename, pathname] = uigetfile('*.mat', 'Load Spectral Hypercube',dir0);
    file_tot=[pathname,filename];

    h=waitbar(1/2,'Loading the Spectra Hypercube');
    input_data=load(file_tot);
    close(h);
    Hypercube=input_data.Hyperspectrum_cube;
    frequency=input_data.fr_real;
    clear input_data.Hyperspectrum_cube input_data.fr_real input_data.f
else
    [Hyperspectrum_cube,f,saturationMap,file_totCal,fr_real,Intens] = H5HypercubeRead;
    Hypercube=Hyperspectrum_cube;
    frequency=fr_real;
    clear Hyperspectrum_cube fr_real f saturationMap file_totCal Intens
end

 n=Dinput('\nNumber of specific frequencies [max 6]:',1);
 
 while n<1 || n>6
      n=Dinput('\nNumber of specific frequencies [max 6]:',1);
 end
 
 pos=zeros(n,1);
 specificFrequencies=zeros(n,1);
 
 fprintf('\nWrite %d specific frequencies in THz',n);
 for k=1:n
     fprintf('\nFrequency [THz] n°%d\n',k);
     specificFrequencies(k)=input('input: ');
     fprintf('\n');
 end
 
 for j=1:n
         target=specificFrequencies(j);
         temp=abs(target-frequency);
         pos(j)=find(temp == min(abs(target - frequency))); %position of the frequency closest to the input one in the frequency array
 end
 
 if ~isempty(pos(pos==0)) %se non vengono trovate le frequenze giuste termina il programma
    fprintf('ERROR. No frequency found.');
    return;
end
 
Intens=sum(abs(Hypercube),3); 
gamma=1;

immagine=zeros(size(Intens,1),size(Intens,2),3);
immagine(:,:,1)=Intens;
immagine(:,:,2)=Intens;
immagine(:,:,3)=Intens;
max_immagine=max(max(max(immagine)));
immagine=immagine./max_immagine;

h0=figure;
image(immagine);
title('Choose gamma correction');

prompt={'Gamma value (-1 and then ENTER to exit)'};
name='Select gamma value';
numlines=[1 35];
defaultanswer={'1'};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

while gamma>0
    prec_gamma=gamma;
    answer=inputdlg(prompt,name,numlines,defaultanswer,options);
    gamma=str2double(cell2mat(answer));

    if gamma>0
        h0;
        image(immagine.^gamma);
        title('Choose gamma correction');
    end
end
pause;
close(h0);
gamma=prec_gamma;
clear prec_gamma;
immagine=immagine.^gamma;

h0=figure;
image(immagine);
title('Zoom the ROI then ENTER');
pause;
V=round(axis);

close(h0);

rmin=round(V(3)); rmax=round(V(4))-1; % First and last row
cmin=round(V(1)); cmax=round(V(2))-1; % First and last column

range_x=rmin:rmax;
range_y=cmin:cmax;
clear cmin cmax

points=zeros(size(range_x,2),size(range_y,2));
meas_frequency=zeros(n,size(points,1)*size(points,2));

for k=1:n
    points=squeeze(Hypercube(range_x,range_y,pos(k)));
    meas_frequency(k,:)=reshape(points,1,size(points,1)*size(points,2));
end

media=zeros(n,1);
r=zeros(n,1);
r_comp=zeros(n,1);
r_diff=zeros(n,1);

for k=1:n
    %calcolo delle medie
    media(k)=mean(squeeze(meas_frequency(k,:)));
    %calcolo delle varianze
    r(k)=sqrt(var(abs(squeeze(meas_frequency(k,:)))));
    r_comp(k)=sqrt(abs(var(squeeze(meas_frequency(k,:)))));
    r_diff(k)=sqrt(var(abs(squeeze(meas_frequency(k,:))-media(k)))); %sigma^2=var|v-<v>|
end
%punti nelle circonferenze
t=linspace(0,2*pi,1000);

figure
MarkerColors = [1.00 0.00 0.00;
             0.00 0.00 1.00;
             0.00 1.00 0.00;
             0.00 1.00 1.00;
             1.00 1.00 0.00;
             1.00 0.00 1.00];
         
%plot dei punti
for k=1:n
    scatter(real(meas_frequency(k,:)),imag(meas_frequency(k,:)),3,'Linewidth',0.1,'MarkerFaceColor',MarkerColors(k,:),'MarkerEdgeColor','w',...
        'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.001);
    hold on
end

circleColor= [0.8 0.0 0.0;
             0.0 0.0 0.8;
             0.0 0.8 0.0;
             0.0 0.8 1.0;
             1.0 0.8 0.0;
             0.8 0.0 1.0];

%plot delle circonferenze
for k=1:n
    plot(real(media(k))+r(k)*cos(t),imag(media(k))+r(k)*sin(t),'Linewidth',2,'Color',circleColor(k,:));
    hold on
    %%% SCOMMENTA SE VUOI I PLOT DEGLI ALTRI TIPI DI ST. DEV.
    % plot(real(media(k))+r_comp(k)*cos(t),imag(media(k))+r_comp(k)*sin(t),'--','Linewidth',2,'Color',circleColor(k,:));
    % hold on
    % plot(real(media(k))+r_diff(k)*cos(t),imag(media(k))+r_diff(k)*sin(t),':','Linewidth',2,'Color',circleColor(k,:));
    % hold on
end

%plot delle medie come vettori
for k=1:n
    quiver(0,0,real(media(k)),imag(media(k)),0,'Linewidth',2,'Color',circleColor(k,:));
end


allValues=zeros(1,size(points,1)*size(points,2)*n);
for k=1:n
    allValues(1,(k-1).*size(points,1)*size(points,2)+1:k.*size(points,1)*size(points,2))=squeeze(meas_frequency(k,:));
end
massimo_real=max(abs(real(allValues)));
massimo_imag=max(abs(imag(allValues)));

axis([-massimo_real-massimo_real/10 massimo_real+massimo_real/10 -massimo_imag-massimo_imag/10 massimo_imag+massimo_imag/10]);
xlabel('Re','FontWeight','bold');
ylabel('Im','FontWeight','bold');
legenda=zeros(n,1);
for k=1:n
    legenda(k)=round(frequency(pos(k)));
end

measUnit=[' THz';' THz';' THz';' THz';' THz';' THz'];
measUnit=measUnit(1:n,:);

legenda=[num2str([legenda]) measUnit];
l=legend({legenda});
set(l,'Position',[0.02 0.78 0.15 0.15]);
grid on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['Gauss Plot, ',num2str(size(range_x,2)),'x',num2str(size(range_y,2)),' selected pixels']);