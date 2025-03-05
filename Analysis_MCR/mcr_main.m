function varargout = mcr_main(varargin)
% MCR_MAIN M-file for mcr_main.fig
%      MCR_MAIN, by itself, creates a new MCR_MAIN or raises the existing
%      singleton*.
%
%      H = MCR_MAIN returns the handle to a new MCR_MAIN or the handle to
%      the existing singleton*.
%
%      MCR_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MCR_MAIN.M with the given input arguments.
%
%      MCR_MAIN('Property','Value',...) creates a new MCR_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mcr_main_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mcr_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mcr_main

% Last Modified by GUIDE v2.5 10-Apr-2015 14:44:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mcr_main_OpeningFcn, ...
    'gui_OutputFcn',  @mcr_main_OutputFcn, ...
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

% --- Executes just before mcr_main is made visible.
function mcr_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mcr_main (see VARARGIN)

% Choose default command line output for mcr_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mcr_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
clc;
disp ('   ')
disp ('   ')
disp ('   ************************************************************************************')
disp ('   *                           MATLAB program MCR-ALS:                                *')
disp ('   *                           Multivariate Curve Resolution (MCR)                    *')
disp ('   *                           Alternating Least Squares (ALS)                        *')
disp ('   ************************************************************************************')
disp ('   ')
disp ('                                                                 MCR 2.0     version   ')
disp ('   ')
disp ('   ')

% initializations
if ((evalin('base','exist(''mcr_als'')'))==1);evalin('base','clear(''mcr_als'')');end;
evalin('base','mcr_als.alsOptions.weighted.constrained=0;');

function varargout = mcr_main_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Data set
% *************************************************************************

% List of the variables of the WS
% *************************************************************************


% --- Executes on button press in check_dataType.
function check_dataType_Callback(hObject, eventdata, handles)
% hObject    handle to check_dataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_dataType

valorDT=get(hObject,'Value');

if valorDT==0
    lsvar=evalin('base','whos');
    [aa,bb]=size(lsvar);
    j=2;
    lsv(1)={'select a variable from the WS'};
    for i=1:aa,
        csize=length(lsvar(i).class);
        if csize==6,
            if lsvar(i).class=='double',
                lsb=[lsvar(i).name];
                lsv(j)={lsb};
                j=j+1;
            end;
        end;
    end;
    
    set(handles.popupmenu3,'string',lsv)
    
else
    lsvar=evalin('base','whos');
    [aa,bb]=size(lsvar);
    j=2;
    lsv(1)={'select a variable from the WS'};
    for i=1:aa,
        csize=length(lsvar(i).class);
        if csize==6,
            if (lsvar(i).class==('double')|lsvar(i).class==('single')|lsvar(i).class==('int8  ')|lsvar(i).class==('int16 ')|lsvar(i).class==('int32 ')|lsvar(i).class==('int64 ')|lsvar(i).class==('uint8 ')|lsvar(i).class==('uint16')|lsvar(i).class==('uint32')|lsvar(i).class==('uint64'));
                lsb=[lsvar(i).name];
                lsv(j)={lsb};
                j=j+1;
            end;
        end;
    end;
    
    set(handles.popupmenu3,'string',lsv)
    
end



% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

lsvar=evalin('base','whos');
[aa,bb]=size(lsvar);
j=2;
lsv(1)={'select a variable from the WS'};
for i=1:aa,
    csize=length(lsvar(i).class);
    if csize==6,
        if lsvar(i).class=='double',
            lsb=[lsvar(i).name];
            lsv(j)={lsb};
            j=j+1;
        end;
    end;
end;

set(hObject,'string',lsv)

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

popmenu3=get(handles.popupmenu3,'String');
pm3=get(handles.popupmenu3,'Value');

if pm3==1
    
else
    selec3=char([popmenu3(pm3)]);
    matdad=evalin('base',selec3);
    assignin('base','matdad',matdad);
    assignin('base','selec3',selec3);
    evalin('base','mcr_als.OriginalData=matdad;');
    evalin('base','mcr_als.data=matdad;');
    evalin('base','mcr_als.aux.data_name=selec3;');
    set(handles.text_nom,'string',selec3);
    evalin('base','clear matdad selec3'); 
    
    % activate: number of components
    set(handles.push_dades,'enable','on');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    evalin('base','mcr_als.aux.estat=0;');
    
    % activate: weight text
    set(handles.text_weight,'string','No');    
end


% Plot the selected data
% *************************************************************************

% --- Executes on button press in push_dades.
function push_dades_Callback(hObject, eventdata, handles)
% hObject    handle to push_dades (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
% 
% figure;
% matriu=evalin('base','mcr_als.data');
% plot(matriu');

plot_dades;


% N components
% *************************************************************************

% --- Executes on button press in but_SVD.
function but_SVD_Callback(hObject, eventdata, handles)
% hObject    handle to but_SVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of but_SVD

set(handles.push_dades,'enable','off');
set(handles.check_weight,'enable','off');
set(handles.push_weight,'enable','off');
set(handles.but_SVD,'enable','off');
set(handles.but_manComp,'enable','off');
set(handles.but_Pure,'enable','off');
set(handles.but_EFA,'enable','off');
set(handles.but_ManEst,'enable','off');
set(handles.but_ALS,'enable','off');
set(handles.but_quit,'enable','off');
evalin('base','mcr_als.aux.estat=1;');

matriu_svd=evalin('base','mcr_als.data');
[a,b]=size(matriu_svd);
evalin('base','mcr_als.aux.estat_svd=1;');

if (a>10000 | b>10000) % avoid problems with big data
    svd_selection;
    while evalin('base','mcr_als.aux.estat_svd;')==1
        pause(0.3);
    end
    nrsvd=evalin('base','mcr_als.aux.nrsvd;');
    [U,S,V]=svds(matriu_svd,nrsvd);
else
    [U,S,V]=svd(matriu_svd,'econ');
end

assignin('base','U',U);
assignin('base','S',S);
assignin('base','V',V);
evalin('base','mcr_als.CompNumb.U=U;');
evalin('base','mcr_als.CompNumb.S=S;');
evalin('base','mcr_als.CompNumb.V=V;');
evalin('base','clear U S V');
svdgui;

while(evalin('base','mcr_als.aux.estat')==1)
    pause(0.3);
end

if (evalin('base','mcr_als.aux.estat')==2)
    set(handles.push_dades,'enable','off');
    set(handles.check_weight,'enable','on');
    set(handles.push_weight,'enable','off');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','on');
    set(handles.but_EFA,'enable','on');
    set(handles.but_ManEst,'enable','on');
    set(handles.but_ALS,'enable','off');
    set(handles.but_quit,'enable','on');
    set(handles.text_nc,'string',num2str(evalin('base','mcr_als.CompNumb.nc')));
    set(handles.text_ncmet,'string',evalin('base','mcr_als.CompNumb.method'));
    
elseif (evalin('base','mcr_als.aux.estat')==3)
    set(handles.push_dades,'enable','off');
    set(handles.check_weight,'enable','off');
    set(handles.push_weight,'enable','off');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','off');
    set(handles.but_EFA,'enable','off');
    set(handles.but_ManEst,'enable','off');
    set(handles.but_ALS,'enable','off');
    set(handles.but_quit,'enable','on');
end

evalin('base','mcr_als.aux.estat=0;');

% --- Executes on button press in but_manComp.
function but_manComp_Callback(hObject, eventdata, handles)
% hObject    handle to but_manComp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of but_manComp

set(handles.push_dades,'enable','off');
set(handles.check_weight,'enable','on');
set(handles.push_weight,'enable','off');
set(handles.but_SVD,'enable','off');
set(handles.but_manComp,'enable','off');
set(handles.but_Pure,'enable','off');
set(handles.but_EFA,'enable','off');
set(handles.but_ManEst,'enable','off');
set(handles.but_ALS,'enable','off');
set(handles.but_quit,'enable','off');
evalin('base','mcr_als.aux.estat=1;');
manComp;

while(evalin('base','mcr_als.aux.estat')==1)
    pause(0.3);
end

if (evalin('base','mcr_als.aux.estat')==2)
    set(handles.push_dades,'enable','off');
    set(handles.check_weight,'enable','on');
    set(handles.push_weight,'enable','off');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','on');
    set(handles.but_EFA,'enable','on');
    set(handles.but_ManEst,'enable','on');
    set(handles.but_ALS,'enable','off');
    set(handles.but_quit,'enable','on');
    set(handles.text_nc,'string',num2str(evalin('base','mcr_als.CompNumb.nc')));
    set(handles.text_ncmet,'string',evalin('base','mcr_als.CompNumb.method'));
elseif (evalin('base','mcr_als.aux.estat')==3)
    set(handles.push_dades,'enable','off');
    set(handles.check_weight,'enable','on');
    set(handles.push_weight,'enable','off');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','off');
    set(handles.but_EFA,'enable','off');
    set(handles.but_ManEst,'enable','off');
    set(handles.but_ALS,'enable','off');
    set(handles.but_quit,'enable','on');
end

evalin('base','mcr_als.aux.estat=0;');


% Weighting
% *************************************************************************

% --- Executes on button press in check_weight.
function check_weight_Callback(hObject, eventdata, handles)
% hObject    handle to check_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_weight

valorW=get(hObject,'Value');

if valorW==0
    set(handles.push_weight,'enable','off');
    set(handles.text_weight,'string','No');
    set(handles.text_ready,'string','-','ForegroundColor','Black');
    set(handles.but_Pure,'enable','on');
    set(handles.but_EFA,'enable','on');
    set(handles.but_ManEst,'enable','on');

    evalin('base','mcr_als.data=mcr_als.OriginalData;'); % reset original data without weighting
    selec3=evalin('base','mcr_als.aux.data_name;');
    set(handles.text_nom,'string',selec3);
else
    set(handles.push_weight,'enable','on');
    set(handles.text_weight,'string','Yes');
    set(handles.text_ready,'string','Not Ready','ForegroundColor','Red');
    set(handles.but_Pure,'enable','off');
    set(handles.but_EFA,'enable','off');
    set(handles.but_ManEst,'enable','off');
end

evalin('base','mcr_als.alsOptions.weighted.appWeight=0;');
assignin('base','valorW',valorW);
evalin('base','mcr_als.alsOptions.weighted.check_application=valorW;');
evalin('base','clear valorW');


% --- Executes on button press in push_weight.
function push_weight_Callback(hObject, eventdata, handles)
% hObject    handle to push_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.push_dades,'enable','inactive');
set(handles.push_weight,'enable','off');
set(handles.but_SVD,'enable','off');
set(handles.but_manComp,'enable','off');
set(handles.but_Pure,'enable','off');
set(handles.but_EFA,'enable','off');
set(handles.but_ManEst,'enable','off');
set(handles.but_ALS,'enable','off');
set(handles.but_quit,'enable','off');

evalin('base','mcr_als.aux.estat=1;');

appWeight=evalin('base','mcr_als.alsOptions.weighted.appWeight;');
if appWeight==0
    evalin('base','mcr_als.alsOptions.weighted.application=0;');
    evalin('base','mcr_als.alsOptions.weighted.max_iterations=200000;');
    evalin('base','mcr_als.alsOptions.weighted.conv_limit=1e-10;');
    evalin('base','mcr_als.alsOptions.weighted.constrained=0;');
    weighted;
else
    weighted;
end

while(evalin('base','mcr_als.aux.estat')==1)
    pause(0.3);
end

if (evalin('base','mcr_als.aux.estat')==2)
    set(handles.push_dades,'enable','off');
    set(handles.push_weight,'enable','on');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','on');
    set(handles.but_EFA,'enable','on');
    set(handles.but_ManEst,'enable','on');
    set(handles.but_ALS,'enable','off');
    set(handles.but_quit,'enable','on');
    if evalin('base','mcr_als.alsOptions.weighted.constrained;')==1;
        set(handles.text_ready,'string','Using & OK','ForegroundColor','Green');
        evalin('base','mcr_als.data=mcr_als.alsOptions.weighted.Xmat_weighted;'); % reset original data without weighting
        selec3=evalin('base','mcr_als.aux.data_name;');
        set(handles.text_nom,'string',[selec3,' - weighted']);
        
    else
        set(handles.text_ready,'string','Done, not OK','ForegroundColor','Blue');
        selec3=evalin('base','mcr_als.aux.data_name;');
        set(handles.text_nom,'string',selec3);
    end
elseif (evalin('base','mcr_als.aux.estat')==3)
    set(handles.push_dades,'enable','off');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','on');
    set(handles.but_EFA,'enable','on');
    set(handles.but_ManEst,'enable','on');
    set(handles.but_ALS,'enable','off');
    set(handles.but_quit,'enable','on');
    set(handles.text_ready,'string','Cancelled','ForegroundColor','Black');
end

evalin('base','mcr_als.aux.estat=0;');


% Initial estimations
% *************************************************************************

% --- Executes on button press in but_Pure.
function but_Pure_Callback(hObject, eventdata, handles)
% hObject    handle to but_Pure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of but_Pure

set(handles.push_dades,'enable','off');
set(handles.push_weight,'enable','off');
set(handles.but_SVD,'enable','off');
set(handles.but_manComp,'enable','off');
set(handles.but_Pure,'enable','off');
set(handles.but_EFA,'enable','off');
set(handles.but_ManEst,'enable','off');
set(handles.but_ALS,'enable','off');
set(handles.but_quit,'enable','off');

evalin('base','mcr_als.aux.estat=1;');

puregui;

while(evalin('base','mcr_als.aux.estat')==1)
    pause(0.3);
end

if (evalin('base','mcr_als.aux.estat')==2)
    set(handles.push_dades,'enable','off');
    set(handles.push_weight,'enable','off');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','on');
    set(handles.but_EFA,'enable','on');
    set(handles.but_ManEst,'enable','on');
    set(handles.but_ALS,'enable','on');
    set(handles.but_quit,'enable','on');
    set(handles.text_iniesta,'string',evalin('base','mcr_als.InitEstim.method'));
elseif (evalin('base','mcr_als.aux.estat')==3)
    set(handles.push_dades,'enable','off');
    set(handles.push_weight,'enable','off');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','on');
    set(handles.but_EFA,'enable','on');
    set(handles.but_ManEst,'enable','on');
    set(handles.but_ALS,'enable','off');
    set(handles.but_quit,'enable','on');
end

evalin('base','mcr_als.aux.estat=0;');

% --- Executes on button press in but_EFA.
function but_EFA_Callback(hObject, eventdata, handles)
% hObject    handle to but_EFA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of but_EFA

set(handles.push_dades,'enable','off');
set(handles.push_weight,'enable','off');
set(handles.but_SVD,'enable','off');
set(handles.but_manComp,'enable','off');
set(handles.but_Pure,'enable','off');
set(handles.but_EFA,'enable','off');
set(handles.but_ManEst,'enable','off');
set(handles.but_ALS,'enable','off');
set(handles.but_quit,'enable','off');

evalin('base','mcr_als.aux.estat=1;');
evalin('base','mcr_als.aux.counter=0;');
efagui;

while(evalin('base','mcr_als.aux.estat')==1)
    pause(0.3);
end

if (evalin('base','mcr_als.aux.estat')==2)
    set(handles.push_dades,'enable','off');
    set(handles.push_weight,'enable','off');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','on');
    set(handles.but_EFA,'enable','on');
    set(handles.but_ManEst,'enable','on');
    set(handles.but_ALS,'enable','on');
    set(handles.but_quit,'enable','on');
    set(handles.text_iniesta,'string',evalin('base','mcr_als.InitEstim.method'));
elseif (evalin('base','mcr_als.aux.estat')==3)
    set(handles.push_dades,'enable','off');
    set(handles.push_weight,'enable','off');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','on');
    set(handles.but_EFA,'enable','on');
    set(handles.but_ManEst,'enable','on');
    set(handles.but_ALS,'enable','off');
    set(handles.but_quit,'enable','on');
end

evalin('base','mcr_als.aux.estat=0;');


% --- Executes on button press in but_ManEst.
function but_ManEst_Callback(hObject, eventdata, handles)
% hObject    handle to but_ManEst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of but_ManEst

set(handles.push_dades,'enable','off');
set(handles.push_weight,'enable','off');
set(handles.but_SVD,'enable','off');
set(handles.but_manComp,'enable','off');
set(handles.but_Pure,'enable','off');
set(handles.but_EFA,'enable','off');
set(handles.but_ManEst,'enable','off');
set(handles.but_ALS,'enable','off');
set(handles.but_quit,'enable','off');

evalin('base','mcr_als.aux.estat=1;');

manEst;

while(evalin('base','mcr_als.aux.estat')==1)
    pause(0.3);
end

if (evalin('base','mcr_als.aux.estat')==2)
    
    set(handles.push_dades,'enable','off');
    set(handles.push_weight,'enable','off');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','on');
    set(handles.but_EFA,'enable','on');
    set(handles.but_ManEst,'enable','on');
    set(handles.but_ALS,'enable','on');
    set(handles.but_quit,'enable','on');
    set(handles.text_iniesta,'string',evalin('base','mcr_als.InitEstim.method'));
elseif (evalin('base','mcr_als.aux.estat')==3)
    
    set(handles.push_dades,'enable','off');
    set(handles.push_weight,'enable','off');
    set(handles.but_SVD,'enable','on');
    set(handles.but_manComp,'enable','on');
    set(handles.but_Pure,'enable','on');
    set(handles.but_EFA,'enable','on');
    set(handles.but_ManEst,'enable','on');
    set(handles.but_ALS,'enable','off');
    set(handles.but_quit,'enable','on');
end

evalin('base','mcr_als.aux.estat=0;');


% Go to ALS
% *************************************************************************

% --- Executes on button press in but_ALS.
function but_ALS_Callback(hObject, eventdata, handles)
% hObject    handle to but_ALS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of but_ALS
close;

% evalin('base','initEstim=mcr_als.InitEstim.iniesta;');

matdad=evalin('base','mcr_als.data;');
mm=max(size(matdad));

if mm>10000;
    answer = questdlg('Warning: large dataset. Graphics?', ...
        'Large dataset warning', ...
        'Yes','No','No');
    % Handle response
    
    switch answer
        case 'Yes'
            als2013;
        case 'No'
            als2013fast;
    end
    
else
    als2013;
end


% Quit
% *************************************************************************

% --- Executes on button press in but_quit.
function but_quit_Callback(hObject, eventdata, handles)
% hObject    handle to but_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of but_quit

if ((evalin('base','exist(''mcr_als'')'))==1);evalin('base','clear(''mcr_als'')');end;
close;



% --- Executes on button press in push_alsInfo.
function push_alsInfo_Callback(hObject, eventdata, handles)
% hObject    handle to push_alsInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

als_info;

