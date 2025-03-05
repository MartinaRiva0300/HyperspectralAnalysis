function varargout = quadrilinear_info(varargin)
% QUADRILINEAR_INFO MATLAB code for quadrilinear_info.fig
%      QUADRILINEAR_INFO, by itself, creates a new QUADRILINEAR_INFO or raises the existing
%      singleton*.
%
%      H = QUADRILINEAR_INFO returns the handle to a new QUADRILINEAR_INFO or the handle to
%      the existing singleton*.
%
%      QUADRILINEAR_INFO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUADRILINEAR_INFO.M with the given input arguments.
%
%      QUADRILINEAR_INFO('Property','Value',...) creates a new QUADRILINEAR_INFO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before quadrilinear_info_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to quadrilinear_info_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help quadrilinear_info

% Last Modified by GUIDE v2.5 08-Jan-2014 13:07:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @quadrilinear_info_OpeningFcn, ...
                   'gui_OutputFcn',  @quadrilinear_info_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before quadrilinear_info is made visible.
function quadrilinear_info_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to quadrilinear_info (see VARARGIN)

% Choose default command line output for quadrilinear_info
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes quadrilinear_info wait for user response (see UIRESUME)
% uiwait(handles.figure1);

    uaug=evalin('base','mcr_als.alsOptions.resultats.optim_concs;');
    vaug=evalin('base','mcr_als.alsOptions.resultats.optim_specs;');
    nmat=evalin('base','mcr_als.alsOptions.trilin.ne;');
    d=evalin('base','mcr_als.alsOptions.matdad;');

% [u1,u2,l,dc,dtc,dt]=umodeaugq(copt,sopt,nmat,1,d);
% umodeaugq

[nrt,ns]=size(uaug);
[ns,nc]=size(vaug);
[nd1,nd2]=size(nmat);
u1=[];u2=[];u3=[];l=[];ll=[];um=[];

if nd2>1,
    disp('quadrilinear model is applied'), 
    iquadril=1;
    %ne1=nmat(1);ne2=nmat(2);ne3=nmat(3);
    ne1=nmat(3);ne2=nmat(2);ne3=nmat(1);
    nr1=nrt/ne1;nr2=nr1/ne2;
    dq=zeros(ne1,ne2,ne3,nc);
    dt=zeros(ne1,ne2*ne3,nc);
end



% *************************************************************************

irow=0;
for i1=1:ne1
    jrow=0;
    for i2=1:ne2
        for i3=1:ne3
            irow=irow+1;
            jrow=jrow+1;
            dq(i1,i2,i3,:)=d(irow,:);
            dt(i1,jrow,:)=d(irow,:);
        end
    end
end

for j=1:ns
    m=[];
    irow=0;
    for i1=1:ne1
        mrow=0;
        for i2=1:ne2*ne3
            irow=irow+1;
            mrow=mrow+1;
            m1(mrow,1)=uaug(irow,j);
        end
        m=[m,m1];
    end
    [u,s,v]=svd(m,0);
    l1=[];for il=1:ne1,l1=[l1;s(il,il)];end
    l=[l,l1];
    % lofbil(m,u(:,1),v(:,1)'.*s(1,1));
    u1=[u1,abs(v(:,1)*s(1,1))];
    um=[um,abs(u(:,1))];
    
    m=[];
    irow=0;
    for i2=1:ne2
        mrow=0;
        for i3=1:ne3
            irow=irow+1;
            mrow=mrow+1;
            m2(mrow,1)=um(irow,j);
        end
        m=[m,m2];
    end

    [u,s,v]=svd(m,0);
    u2=[u2,abs(v(:,1)*s(1,1))];
    u3=[u3,abs(u(:,1))];
    ll1=[];for il=1:ne2,ll1=[ll1;s(il,il)];end
    ll=[ll,ll1];
    size(u2);
    size(u3);
    %pause
        
% repeat for all components, ns

end


disp('check of multilinearity from svd of u1 folded profiles: ');disp(l)
disp('check of multilinearity from svd of u2-u3 folded profiles: ');disp(ll)

% Test lofbil (x,uaug,vaug) 

disp(' ');disp('lof and R2 for a bilinear model')
ifn=find(isfinite(d)==1);
for i=1:nrt
    for j=1:nc
        if isfinite(d(i,j))==1,
            dc(i,j)=0;
            for ls=1:ns
                dc(i,j)=dc(i,j)+uaug(i,ls)*vaug(ls,j);
            end
        end
    end
end

res=d(ifn)-dc(ifn);
sumd=sum(sum(d(ifn).*d(ifn)));
sumdc=sum(sum(dc(ifn).*dc(ifn)));
sumres=sum(sum(res.*res));
% disp(['sstot,sscalc and ssres = ',num2str([sumd,sumdc,sumres])]);
loft=(sqrt(sumres/sumd))*100;
r2t1=((sumdc)/sumd)*100;
r2t2=((sumd-sumres)/sumd)*100;
disp(['lof (%) = ',num2str(loft)]);
disp(['R2 % = ',num2str(r2t2)]);

set(handles.text_2LOF,'string',num2str(loft));
set(handles.text_2R2,'string',num2str(r2t2));

% Test loftril (x,u1,um=u2*u3,vaug')

disp(' ');disp('lof and R2 for a trilinear model')
for i=1:ne1
        for j=1:ne2*ne3
            for k=1:nc
                if isfinite(dt(i,j,k))==1,
                    dtc(i,j,k)=0;
                    for ls=1:ns
                        dtc(i,j,k)=dtc(i,j,k)+u1(i,ls)*um(j,ls)*vaug(ls,k);
                    end
                end
            end
        end
end


ifn=find(isfinite(dt)==1);
res=dt(ifn)-dtc(ifn);
sumdt=sum(sum(sum(dt(ifn).*dt(ifn))));
sumdtc=sum(sum(sum(dtc(ifn).*dtc(ifn))));
sumres=sum(sum(sum(res.*res)));
% disp(['sstot,sscalc and ssres = ',num2str([sumdt,sumdtc,sumres])]);
loft=(sqrt(sumres/sumdt))*100;
r2t1=((sumdtc)/sumdt)*100;
r2t2=((sumdt-sumres)/sumdt)*100;
disp(['lof (%) = ',num2str(loft)]);
disp(['R2 % = ',num2str(r2t2)]);

set(handles.text_3LOF,'string',num2str(loft));
set(handles.text_3R2,'string',num2str(r2t2));


% Test lofquadril (x,u1,u2,u3,vaug')

disp(' ');disp('lof and R2 for a quadrilinear model')
for i=1:ne1
    for j=1:ne2
        for k=1:ne3
            for l=1:nc
                if isfinite(dq(i,j,k,l))==1,
                    dqc(i,j,k,l)=0;
                    for ls=1:ns
                        dqc(i,j,k,l)=dqc(i,j,k,l)+u1(i,ls)*u2(j,ls)*u3(k,ls)*vaug(ls,l);
                    end
                end
            end
        end
    end
end

% Build the column-wise augmented recalculated two-way data matrix

nrow=0;
for i=1:ne1
    for j=1:ne2
        for k=1:ne3
            %for l=1:ne4
                nrow=nrow+1;
                for ls=1:ns
                        dc(nrow,:)=dqc(i,j,k,:);
                end
            % end
        end
    end
end

ifn=find(isfinite(dq)==1);
res=dq(ifn)-dqc(ifn);
sumdt=sum(sum(sum(sum(dq(ifn).*dq(ifn)))));
sumdtc=sum(sum(sum(sum(dqc(ifn).*dqc(ifn)))));
sumres=sum(sum(sum(sum((res.*res)))));
% disp(['sstot,sscalc and ssres = ',num2str([sumdt,sumdtc,sumres])]);
loft=(sqrt(sumres/sumdt))*100;
r2t1=((sumdtc)/sumdt)*100;
r2t2=((sumdt-sumres)/sumdt)*100;
disp(['lof (%) = ',num2str(loft)]);
disp(['R2 % = ',num2str(r2t2)]);

set(handles.text_4LOF,'string',num2str(loft));
set(handles.text_4R2,'string',num2str(r2t2));



%save umodeaug variables
% [u1,u2,l,dc,dtc,dt]=umodeaug(copt,sopt,nmat,1,);

assignin('base','u1QUAT',u1);
evalin('base','mcr_als.aux.u1=u1QUAT;');
evalin('base','clear u1QUAT');

assignin('base','u2QUAT',u2);
evalin('base','mcr_als.aux.u2=u2QUAT;');
evalin('base','clear u2QUAT');

assignin('base','u3QUAT',u3);
evalin('base','mcr_als.aux.u3=u3QUAT;');
evalin('base','clear u3QUAT');

assignin('base','vaugQUAT',vaug);
evalin('base','mcr_als.aux.vaug=vaugQUAT;');
evalin('base','clear vaugQUAT');

assignin('base','lQUAT',l);
evalin('base','mcr_als.aux.l=lQUAT;');
evalin('base','clear lQUAT');

assignin('base','llQUAT',ll);
evalin('base','mcr_als.aux.ll=llQUAT;');
evalin('base','clear llQUAT');

assignin('base','dcQUAT',dc);
evalin('base','mcr_als.aux.dc=dcQUAT;');
evalin('base','clear dcQUAT');

assignin('base','dqcQUAT',dqc);
evalin('base','mcr_als.aux.dqc=dqcQUAT;');
evalin('base','clear dqcQUAT');

assignin('base','dqTRIL',dq);
evalin('base','mcr_als.aux.dq=dqTRIL;');
evalin('base','clear dqTRIL');



% --- Outputs from this function are returned to the command line.
function varargout = quadrilinear_info_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function popup_species_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

nc=evalin('base','mcr_als.CompNumb.nc');
j=2;
llistannc(1)={'select a species...'};
for i=1:nc;
    llistnnc=['Species nr. ',num2str(i)];
    llistannc(j)={llistnnc};
    j=j+1;
end;
set(hObject,'string',llistannc)

% --- Executes on selection change in popup_species.
function popup_species_Callback(hObject, eventdata, handles)
% hObject    handle to popup_species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_species contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_species

valorpop=get(hObject,'Value')-1;

if valorpop>0
    
    u1=evalin('base','mcr_als.aux.u1;');
    u2=evalin('base','mcr_als.aux.u2;');
    u3=evalin('base','mcr_als.aux.u3;');   
    vaug=evalin('base','mcr_als.aux.vaug;');
    
    % mode 1
    axes(handles.axes4);
    plot(u1(:,valorpop));
    mm=length(u1);
    axis([0 mm+1 min(u1(:,valorpop))-0.1*min(u1(:,valorpop)) max(u1(:,valorpop))+0.1*max(u1(:,valorpop))]);
    %title('Mode 1');
    title('Mode 4');

    % mode 2
    axes(handles.axes3);
    plot(u2(:,valorpop));
    mm=length(u2);
    axis([0 mm+1 min(u2(:,valorpop))-0.1*min(u2(:,valorpop)) max(u2(:,valorpop))+0.1*max(u2(:,valorpop))]);
    % title('Mode 2');  
    title('Mode 3');  
    
    % mode 1
    axes(handles.axes1);
    plot(u3(:,valorpop));
    mm=length(u3);
    axis([0 mm+1 min(u3(:,valorpop))-0.1*min(u3(:,valorpop)) max(u3(:,valorpop))+0.1*max(u3(:,valorpop))]);
    % title('Mode 3');
    title('Mode 1');
    
    % mode 4
    axes(handles.axes2);
    plot(vaug(valorpop,:));
    mm=length(vaug);
    axis([0 mm+1 min(vaug(valorpop,:))-0.1*min(vaug(valorpop,:)) max(vaug(valorpop,:))+0.1*max(vaug(valorpop,:))]);
    %title('Mode 4');
    title('Mode 2');
else
    
    axes(handles.axes1);
    cla;
    axes(handles.axes2);
    cla;
    axes(handles.axes3);
    cla;
    axes(handles.axes4);
    cla;
end

% --- Executes on button press in push_export.
function push_export_Callback(hObject, eventdata, handles)
% hObject    handle to push_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

u1=evalin('base','mcr_als.aux.u1;');
assignin('base','u1',u1);

u2=evalin('base','mcr_als.aux.u2;');
assignin('base','u2',u2);

u3=evalin('base','mcr_als.aux.u3;');
assignin('base','u3',u3);

vaug=evalin('base','mcr_als.aux.vaug;');
assignin('base','vaug',vaug);

% --- Executes on button press in push_multiplot.
function push_multiplot_Callback(hObject, eventdata, handles)
% hObject    handle to push_multiplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

kk=0;
ns=evalin('base','mcr_als.CompNumb.nc');

u1=evalin('base','mcr_als.aux.u1;');
u2=evalin('base','mcr_als.aux.u2;');
u3=evalin('base','mcr_als.aux.u3;');
vaug=evalin('base','mcr_als.alsOptions.resultats.optim_specs;');

figure('Name','quadrilinear info multiplot');

for k=1:ns
    if k==1
        kk=kk+1;subplot(ns,4,kk),plot(u3(:,k),'c'),title('Mode 1');ylabel('Component 1');
        kk=kk+1;subplot(ns,4,kk),plot(vaug(k,:),'b'),title('Mode 2');
        kk=kk+1;subplot(ns,4,kk),plot(u2(:,k),'g'),title('Mode 3');
        kk=kk+1;subplot(ns,4,kk),plot(u1(:,k),'r'),title('Mode 4');
    else
            kk=kk+1;subplot(ns,4,kk),plot(u3(:,k),'c');ylabel(['Component ',num2str(k)]);
            kk=kk+1;subplot(ns,4,kk),plot(vaug(k,:),'b')
            kk=kk+1;subplot(ns,4,kk),plot(u2(:,k),'g')
            kk=kk+1;subplot(ns,4,kk),plot(u1(:,k),'r')
    end
end



% --- Executes on button press in push_close.
function push_close_Callback(hObject, eventdata, handles)
% hObject    handle to push_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close;
