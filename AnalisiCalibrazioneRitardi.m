% Per la 633

% Caricamento del dato
% uiopen('/home/cristian/Experiments/Wedges/2022/DelayCharacterization/mean HyperMatrix_633.fig',1)
openfig('C:/Users/HARDi/Desktop/Tesi/Codice immagini/220504/633_19um_0/mean HyperMatrix_633.fig');
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
 
x633=(get(dataObjs,'XData'));
y633=(get(dataObjs,'YData'));

% Per la 633
fr=linspace(-0.015,0.015,20001);
y633_1=y633.*exp(1j*2*pi*0.0611745.*x633); % Demodulo alla frequenza del modo
Y1=FourierDir(x633,y633_1,fr);
plot(fr,abs(Y1))

y633_2=FourierInv(fr,Y1,x633);
plot(x633,unwrap(angle(y633_2)));

ph633=unwrap(angle(y633_2));
dx633=ph633/(2*pi*0.0611745);

fr=linspace(0.055,0.065,1000);
apo=Apodization(1,3599,3599/2);
YF633_def=FourierDir(x633-dx633,y633.*apo,fr);
YF633=FourierDir(x633,y633.*apo,fr);

figure(10);
plot(fr,abs(YF633),'r',fr,abs(YF633_def),'b','linewidth',3)
hold on;

% Per la 543
% Caricamento del dato
% uiopen('/home/cristian/Experiments/Wedges/2022/DelayCharacterization/mean HyperMatrix_543.fig',1)
dir='C:\Users\HARDi\Desktop\Tesi\Codice immagini\220504\543_19um_0\';
openfig([dir 'meanHyperMatrix_543_19um_0.fig']);
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
 
x543=(get(dataObjs,'XData'));
y543=(get(dataObjs,'YData'));
n_ele=length(x543);

y543=y543-mean(y543);

fr=linspace(0.031,0.1,20001);
% fr=linspace(0.018,0.025,2001);
Y1=FourierDir(x543,y543,fr);

% Analisi del dato

fr=linspace(-0.001,0.001,20001);
y543_1=y543.*exp(1j*2*pi*0.073934.*x543); % Demodulo alla frequenza esatta del modo % era 0.074203 % era 0.073939 0.0213075
Y1=FourierDir(x543,y543_1,fr);
% plot(fr,abs(Y1))

y543_2=FourierInv(fr,Y1,x543);
% plot(x543,unwrap(angle(y543_2)));

ph543=unwrap(angle(y543_2));
dx543=ph543/(2*pi*0.073934); % Era 0.074203 % 0.073939 % 0.0213075

% fr=linspace(0.018,0.025,1000);
fr=linspace(0.059,0.085,1000);
apo=Apodization(1,n_ele,n_ele/2);
YF543_def=FourierDir(x543-dx543,y543.*apo,fr); % era -
YF543=FourierDir(x543,y543.*apo,fr);

figure(1)
plot(fr,abs(YF543),'r',fr,abs(YF543_def),'b','linewidth',3)

dxmean=(dx633+dx543)/2;
dxmean=dxmean-mean(dxmean);

fr=linspace(0.055,0.065,1000);
apo=Apodization(1,3599,3599/2);
YF633_def1=FourierDir(x633-dxmean,y633.*apo,fr);
YF633=FourierDir(x633,y633.*apo,fr);

figure(10);
plot(fr,abs(YF633),'r',fr,abs(YF633_def1),'b','linewidth',3);
hold on

fr=linspace(0.065,0.085,1000);
apo=Apodization(1,3599,3599/2);
YF543_def1=FourierDir(x543-dxmean,y543.*apo,fr);
YF543=FourierDir(x543,y543.*apo,fr);

plot(fr,abs(YF543),'r',fr,abs(YF543_def1),'b','linewidth',3)