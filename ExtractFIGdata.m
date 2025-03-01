% Extracts data from a 2D surface plot

h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes

DD=dataObjs{2};
xx0=cell2mat(get(DD,'XData'));
yy0=cell2mat(get(DD,'YData'));

% oppure

xx1=cell2mat(get(dataObjs,'XData'));
yy1=cell2mat(get(dataObjs,'YData'));

xx1=(get(dataObjs,'XData'));
yy1=(get(dataObjs,'YData'));


% oppure
xx0=(get(DD,'XData'));
yy0=(get(DD,'YData'));
zz0=(get(DD,'CData'));