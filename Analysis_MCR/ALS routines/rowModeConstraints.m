function varargout = rowModeConstraints(varargin)
% ROWMODECONSTRAINTS MATLAB code for rowModeConstraints.fig
%      ROWMODECONSTRAINTS, by itself, creates a new ROWMODECONSTRAINTS or raises the existing
%      singleton*.
%
%      H = ROWMODECONSTRAINTS returns the handle to a new ROWMODECONSTRAINTS or the handle to
%      the existing singleton*.
%
%      ROWMODECONSTRAINTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROWMODECONSTRAINTS.M with the given input arguments.
%
%      ROWMODECONSTRAINTS('Property','Value',...) creates a new ROWMODECONSTRAINTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rowModeConstraints_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rowModeConstraints_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rowModeConstraints

% Last Modified by GUIDE v2.5 20-Nov-2012 21:46:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rowModeConstraints_OpeningFcn, ...
                   'gui_OutputFcn',  @rowModeConstraints_OutputFcn, ...
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

% --- Executes just before rowModeConstraints is made visible.
function rowModeConstraints_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rowModeConstraints (see VARARGIN)

% Choose default command line output for rowModeConstraints
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rowModeConstraints wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = rowModeConstraints_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% Non-negativity
% *********************************************************************

% --- Executes on button press in check_nnC.
function check_nnC_Callback(hObject, eventdata, handles)
% hObject    handle to check_nnC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_nnC

noneg=get(hObject,'Value');

if noneg==0;
    set(handles.text_nnCi,'enable','off');
    set(handles.popup_nnCi,'enable','off','value',1);
    set(handles.text_nnCsp,'enable','off'); 
    set(handles.popup_nnCsp,'enable','off','value',1);
    set(handles.text_nnCp,'enable','off');
    set(handles.edit_nnCp,'enable','off','string',' ');
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''nonegC'');');
else
    set(handles.text_nnCi,'enable','on');
    set(handles.popup_nnCi,'enable','on','value',1);
    set(handles.text_nnCsp,'enable','off'); 
    set(handles.popup_nnCsp,'enable','off','value',1);
    set(handles.text_nnCp,'enable','off');
    set(handles.edit_nnCp,'enable','off','string',' ');
end;

assignin('base','noneg',noneg);
evalin('base','mcr_als.alsOptions.nonegC.noneg=noneg;');
evalin('base','clear noneg');

% --- Executes during object creation, after setting all properties.
function popup_nnCi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_nnCi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

llista(1)={'select...'};
llista(2)={'forced to zero'};
llista(3)={'nnls'};
llista(4)={'fnnls'};
set(hObject,'string',llista)

% --- Executes on selection change in popup_nnCi.
function popup_nnCi_Callback(hObject, eventdata, handles)
% hObject    handle to popup_nnCi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_nnCi contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_nnCi

ialg=get(hObject,'Value')-2; 
assignin('base','ialg',ialg);
evalin('base','mcr_als.alsOptions.nonegC.ialg=ialg;');

if ialg==-1
    
    set(handles.text_nnCsp,'enable','off');
    set(handles.popup_nnCsp,'enable','off'),
    set(handles.popup_nnCsp,'string','select...');
    
    set(handles.text_nnCp,'enable','off');
    set(handles.edit_nnCp,'enable','off','string','');
    
elseif ialg==0
    
    dim=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    j=2;
    llistannc(1)={'select...'};
    for i=0:1:dim;
        llistnnc=[i];
        llistannc(j)={llistnnc};
        j=j+1;
    end;
    
    set(handles.text_nnCsp,'enable','on');
    set(handles.popup_nnCsp,'enable','on','value',1)
    set(handles.popup_nnCsp,'string',llistannc)
    
    set(handles.text_nnCp,'enable','off');
    set(handles.edit_nnCp,'enable','off','string','');
       
elseif (ialg==1 | ialg==2)
    dim=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    llistannc(1)={'select...'};
    llistnnc=[dim];
    llistannc(2)={llistnnc};
    
    set(handles.text_nnCsp,'enable','on');
    set(handles.popup_nnCsp,'enable','on','value',1)
    set(handles.popup_nnCsp,'string',llistannc)
    
    set(handles.text_nnCp,'enable','off');
    set(handles.edit_nnCp,'enable','off','string',''); 
    
end

evalin('base','clear ialg');


% --- Executes during object creation, after setting all properties.
function popup_nnCsp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_nnCsp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popup_nnCsp.
function popup_nnCsp_Callback(hObject, eventdata, handles)
% hObject    handle to popup_nnCsp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_nnCsp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_nnCsp

ncneg=get(hObject,'Value')-2;   

assignin('base','ncneg',ncneg);
evalin('base','mcr_als.alsOptions.nonegC.ncneg=ncneg;');

dim=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
nexp=evalin('base','mcr_als.alsOptions.nexp');
ialg=evalin('base','mcr_als.alsOptions.nonegC.ialg');

if ialg == 0
    if ncneg == dim;
        cneg = ones(nexp,dim);
        set(handles.text_nnCp,'enable','off');
        set(handles.edit_nnCp,'enable','off','string',''); 
    elseif  ncneg == -1
        cneg=zeros(nexp,dim);
        set(handles.text_nnCp,'enable','off');
        set(handles.edit_nnCp,'enable','off','string','');
    elseif  ncneg == 0
        cneg = zeros(nexp,dim);
        set(handles.text_nnCp,'enable','off');
        set(handles.edit_nnCp,'enable','off','string','');   
    else
        cneg=[];
        set(handles.text_nnCp,'enable','on');
        set(handles.edit_nnCp,'enable','on'); 
    end
else
%     cneg = ones(nexp,dim);
    if ncneg<0;
        form_cneg = zeros(1,dim);
        cneg(1,:)=form_cneg;
    else;
        form_cneg= ones(1,dim);
        cneg(1,:)=form_cneg;
    end
        
end

assignin('base','cneg',cneg);
evalin('base','mcr_als.alsOptions.nonegC.cneg=cneg;');
evalin('base','clear cneg ncneg');


% --- Executes during object creation, after setting all properties.
function edit_nnCp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nnCp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_nnCp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nnCp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nnCp as text
%        str2double(get(hObject,'String')) returns contents of edit_nnCp as a double

cnegr=str2num(get(hObject,'String'));
cneg = ones(evalin('base','mcr_als.alsOptions.nexp'),1)*cnegr;
assignin('base','cneg',cneg);
evalin('base','mcr_als.alsOptions.nonegC.cneg=cneg;');
evalin('base','clear cneg');


% Unimodality
% *********************************************************************

% --- Executes on button press in check_unimod.
function check_unimod_Callback(hObject, eventdata, handles)
% hObject    handle to check_unimod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_unimod


unimodal=get(hObject,'Value');

if unimodal==1;
    set(handles.text_Ui,'enable','on');
    set(handles.popup_Ui,'enable','on','value',1);     
    set(handles.text_Up,'enable','on');
    set(handles.popup_Up,'enable','on','value',1);    
    set(handles.text_Utol,'enable','on');
    set(handles.edit_Utol,'enable','on','string',' 1.1');   
    set(handles.text_Uv,'enable','off');
    set(handles.edit_Uv,'enable','off','string',' ');
    
    % default value of Utol
    rmod=1.1;
    assignin('base','rmod',rmod);
    evalin('base','mcr_als.alsOptions.unimodC.rmod=rmod;');
    evalin('base','clear rmod');
    
 else
    set(handles.text_Ui,'enable','off');
    set(handles.popup_Ui,'enable','off','value',1);   
    set(handles.text_Up,'enable','off');
    set(handles.popup_Up,'enable','off','value',1);    
    set(handles.text_Utol,'enable','off');
    set(handles.edit_Utol,'enable','off','string',' ');   
    set(handles.text_Uv,'enable','off');
    set(handles.edit_Uv,'enable','off','string',' ');   
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''unimodC'');');  
end;

assignin('base','unimodal',unimodal);
evalin('base','mcr_als.alsOptions.unimodC.unimodal=unimodal;');
evalin('base','clear unimodal');



% --- Executes during object creation, after setting all properties.
function popup_Ui_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Ui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

llista(1)={'select...'};
llista(2)={'vertical'};
llista(3)={'horizontal'};
llista(4)={'average'};
set(hObject,'string',llista)

% --- Executes on selection change in popup_Ui.
function popup_Ui_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Ui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Ui contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Ui

cmod=get(hObject,'Value')-2; 
assignin('base','cmod',cmod);
evalin('base','mcr_als.alsOptions.unimodC.cmod=cmod;');
evalin('base','clear cmod');

% --- Executes during object creation, after setting all properties.
function popup_Up_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

dim=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
j=2;
llistannc(1)={'select...'};
for i=0:1:dim;
    llistnnc=[i];
    llistannc(j)={llistnnc};
    j=j+1;
end;

set(hObject,'string',llistannc)


% --- Executes on selection change in popup_Up.
function popup_Up_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Up contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Up

nmod=get(hObject,'Value')-2; 

dim=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
nexp=evalin('base','mcr_als.alsOptions.nexp');

if nmod == dim
    spmod=ones(1,dim);
    set(handles.text_Uv,'enable','off');
    set(handles.edit_Uv,'enable','off','string','');
elseif  nmod == 0
    spmod=zeros(1,dim);
    set(handles.text_Uv,'enable','off');
    set(handles.edit_Uv,'enable','off','string','');   
elseif  nmod == -1
    spmod=[];
    set(handles.text_Uv,'enable','off');
    set(handles.edit_Uv,'enable','off','string','');
else
    spmod=[];
    set(handles.text_Uv,'enable','on');
    set(handles.edit_Uv,'enable','on'); 
end

assignin('base','spmod',spmod);
evalin('base','mcr_als.alsOptions.unimodC.spmod=spmod;');
evalin('base','clear spmod');


% --- Executes during object creation, after setting all properties.
function edit_Uv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Uv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Uv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Uv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Uv as text
%        str2double(get(hObject,'String')) returns contents of edit_Uv as a double

spmod= str2num(get(hObject,'String'));
assignin('base','spmod',spmod);
evalin('base','mcr_als.alsOptions.unimodC.spmod=spmod;');
evalin('base','clear spmod');


% --- Executes during object creation, after setting all properties.
function edit_Utol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Utol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Utol_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Utol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Utol as text
%        str2double(get(hObject,'String')) returns contents of edit_Utol as a double

rmod= str2num(get(hObject,'String'));

if rmod==1,
    rmod=1.0001;
end

assignin('base','rmod',rmod);
evalin('base','mcr_als.alsOptions.unimodC.rmod=rmod;');
evalin('base','clear rmod');


% Closure
% *********************************************************************

% --- Executes on button press in check_closC.
function check_closC_Callback(hObject, eventdata, handles)
% hObject    handle to check_closC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_closC

closure=get(hObject,'Value');

if closure==1;
    
    set(handles.text_Cnr,'enable','on');
    set(handles.popup_Cnr,'enable','on','value',1);
    
    assignin('base','dc',1);
    evalin('base','mcr_als.alsOptions.closure.dc=dc;');
    evalin('base','clear dc');
    
else
    
    set(handles.text_Cnr,'enable','off');
    set(handles.popup_Cnr,'enable','off','value',1); 
    set(handles.check_Cvar,'enable','off','value',0);
    
    set(handles.text_C1,'enable','off');
    set(handles.text_C1v,'enable','off');
    set(handles.text_C1c,'enable','off');
    set(handles.text_C1sp,'enable','off');
    set(handles.check_C1sp,'enable','off','value',0);
    set(handles.popup_C1c,'enable','off','value',1);
    set(handles.edit_C1,'enable','off','string',' ');
    set(handles.edit_C1v,'enable','off','string',' ');
    set(handles.edit_C1sp,'enable','off','string',' ');
    
    set(handles.text_C2,'enable','off');
    set(handles.text_C2v,'enable','off');
    set(handles.text_C2c,'enable','off');
    set(handles.text_C2sp,'enable','off');
    set(handles.check_C2sp,'enable','off','value',0);
    set(handles.popup_C2c,'enable','off','value',1);
    set(handles.edit_C2,'enable','off','string',' ');
    set(handles.edit_C2v,'enable','off','string',' ');
    set(handles.edit_C2sp,'enable','off','string',' ');    
       
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''closure'');');
    closX='None';
    assignin('base','closX',closX);
    evalin('base','mcr_als.alsOptions.closure.type=closX;');
    evalin('base','clear closX');
    
end;

assignin('base','closure',closure);
evalin('base','mcr_als.alsOptions.closure.closure=closure;');
evalin('base','clear closure');

vc=0;
assignin('base','vc',vc);
evalin('base','mcr_als.alsOptions.closure.vc=vc;');
evalin('base','clear vc');

% --- Executes during object creation, after setting all properties.
function popup_Cnr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Cnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

llista(1)={'select...'};
llista(2)={'1'};
llista(3)={'2'};
set(hObject,'string',llista)

% --- Executes on selection change in popup_Cnr.
function popup_Cnr_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Cnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Cnr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Cnr

iclos=get(hObject,'Value')-1; 

if iclos == 1;
    set(handles.check_Cvar,'enable','on','value',0);
    
    set(handles.text_C1,'enable','on');
    set(handles.edit_C1,'enable','on');
    set(handles.text_C1c,'enable','on');
    set(handles.popup_C1c,'enable','on','value',1);
    set(handles.text_C1sp,'enable','on');
    set(handles.edit_C1sp,'enable','on');
    set(handles.check_C1sp,'enable','on','value',0);
    set(handles.text_C1v,'enable','off');
    set(handles.edit_C1v,'enable','off','string',' ');
    
    set(handles.text_C2,'enable','off');
    set(handles.edit_C2,'enable','off','string',' ');
    set(handles.text_C2c,'enable','off');
    set(handles.popup_C2c,'enable','off','value',1);
    set(handles.text_C2sp,'enable','off');
    set(handles.edit_C2sp,'enable','off','string',' ');
    set(handles.check_C2sp,'enable','off','value',0);
    set(handles.text_C2v,'enable','off');
    set(handles.edit_C2v,'enable','off','string',' ');
    
elseif  iclos==2
    warndlg('Warning: there should not be common species to the two closures','WARNING!!');
    
    set(handles.check_Cvar,'enable','on','value',0);
    
    set(handles.text_C1,'enable','on');
    set(handles.edit_C1,'enable','on');
    set(handles.text_C1c,'enable','on');
    set(handles.popup_C1c,'enable','on','value',1);
    set(handles.text_C1sp,'enable','on');
    set(handles.edit_C1sp,'enable','on');
    set(handles.check_C1sp,'enable','on','value',0);
    set(handles.text_C1v,'enable','off');
    set(handles.edit_C1v,'enable','off','string',' ');
    
    set(handles.text_C2,'enable','on');
    set(handles.edit_C2,'enable','on');
    set(handles.text_C2c,'enable','on');
    set(handles.popup_C2c,'enable','on','value',1);
    set(handles.text_C2sp,'enable','on');
    set(handles.edit_C2sp,'enable','on');
    set(handles.check_C2sp,'enable','on','value',0);
    set(handles.text_C2v,'enable','off');
    set(handles.edit_C2v,'enable','off','string',' ');
else
    set(handles.check_Cvar,'enable','off','value',0);
    
    set(handles.text_C1,'enable','off');
    set(handles.edit_C1,'enable','off');
    set(handles.text_C1c,'enable','off');
    set(handles.popup_C1c,'enable','off','value',1);
    set(handles.text_C1sp,'enable','off');
    set(handles.edit_C1sp,'enable','off');
    set(handles.check_C1sp,'enable','off','value',0);
    set(handles.text_C1v,'enable','off');
    set(handles.edit_C1v,'enable','off','string',' ');
    
    set(handles.text_C2,'enable','off');
    set(handles.edit_C2,'enable','off');
    set(handles.text_C2c,'enable','off');
    set(handles.popup_C2c,'enable','off','value',1);
    set(handles.text_C2sp,'enable','off');
    set(handles.edit_C2sp,'enable','off');
    set(handles.check_C2sp,'enable','off','value',0);
    set(handles.text_C2v,'enable','off');
    set(handles.edit_C2v,'enable','off','string',' ');
end

assignin('base','iclos',iclos);
evalin('base','mcr_als.alsOptions.closure.iclos=iclos;');
evalin('base','clear iclos');

% --- Executes during object creation, after setting all properties.
function edit_C1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_C1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_C1 as text
%        str2double(get(hObject,'String')) returns contents of edit_C1 as a double

tclos1= str2num(get(hObject,'String'));
assignin('base','tclos1',tclos1);
evalin('base','mcr_als.alsOptions.closure.tclos1=tclos1;');
evalin('base','clear tclos1');

% --- Executes during object creation, after setting all properties.
function popup_C1c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_C1c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

llista(1)={'select...'};
llista(2)={'equal to'};
llista(3)={'least-squares closure'};
llista(4)={'lower than or equal to'};
set(hObject,'string',llista)

% --- Executes on selection change in popup_C1c.
function popup_C1c_Callback(hObject, eventdata, handles)
% hObject    handle to popup_C1c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_C1c contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_C1c

iclos1=get(hObject,'Value')-1; 
assignin('base','iclos1',iclos1);
evalin('base','mcr_als.alsOptions.closure.iclos1=iclos1;');
evalin('base','clear iclos1');

% --- Executes on button press in check_C1sp.
function check_C1sp_Callback(hObject, eventdata, handles)
% hObject    handle to check_C1sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_C1sp

check1clos=get(hObject,'Value');

if check1clos==0
    nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    ceros=zeros(1,nsign);
    set(handles.text_C1sp,'enable','on');
    set(handles.edit_C1sp,'enable','on','string',' ');
    assignin('base','sclos1',ceros);
    evalin('base','mcr_als.alsOptions.closure.sclos1=sclos1;');
    evalin('base','clear sclos1');
    
elseif check1clos==1
    nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    unos=ones(1,nsign);    
    set(handles.text_C1sp,'enable','on');
    set(handles.edit_C1sp,'enable','inactive','string',num2str(unos));
    assignin('base','sclos1',unos);
    evalin('base','mcr_als.alsOptions.closure.sclos1=sclos1;');
    evalin('base','clear sclos1');
end



% --- Executes during object creation, after setting all properties.
function edit_C1sp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_C1sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_C1sp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_C1sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_C1sp as text
%        str2double(get(hObject,'String')) returns contents of edit_C1sp as a double

sclos1= str2num(get(hObject,'String'));
assignin('base','sclos1',sclos1);
evalin('base','mcr_als.alsOptions.closure.sclos1=sclos1;');
evalin('base','clear sclos1');

% --- Executes during object creation, after setting all properties.
function edit_C2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_C2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_C2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_C2 as text
%        str2double(get(hObject,'String')) returns contents of edit_C2 as a double

tclos2= str2num(get(hObject,'String'));
assignin('base','tclos2',tclos2);
evalin('base','mcr_als.alsOptions.closure.tclos2=tclos2;');
evalin('base','clear tclos2');


% --- Executes during object creation, after setting all properties.
function popup_C2c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_C2c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

llista(1)={'select...'};
llista(2)={'equal to'};
llista(3)={'least-squares closure'};
llista(4)={'lower than or equal to'};
set(hObject,'string',llista)

% --- Executes on selection change in popup_C2c.
function popup_C2c_Callback(hObject, eventdata, handles)
% hObject    handle to popup_C2c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_C2c contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_C2c

iclos2=get(hObject,'Value')-1; 
assignin('base','iclos2',iclos2);
evalin('base','mcr_als.alsOptions.closure.iclos2=iclos2;');
evalin('base','clear iclos2');

% --- Executes on button press in check_C2sp.
function check_C2sp_Callback(hObject, eventdata, handles)
% hObject    handle to check_C2sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_C2sp

check2clos=get(hObject,'Value');

if check2clos==0
    nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    ceros=zeros(1,nsign);
    set(handles.text_C2sp,'enable','on');
    set(handles.edit_C2sp,'enable','on','string',' ');
    assignin('base','sclos2',ceros);
    evalin('base','mcr_als.alsOptions.closure.sclos2=sclos2;');
    evalin('base','clear sclos2');
    
elseif check2clos==1
    nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    unos=ones(1,nsign);    
    set(handles.text_C2sp,'enable','on');
    set(handles.edit_C2sp,'enable','inactive','string',num2str(unos));
    assignin('base','sclos2',unos);
    evalin('base','mcr_als.alsOptions.closure.sclos2=sclos2;');
    evalin('base','clear sclos2');
end


% --- Executes during object creation, after setting all properties.
function edit_C2sp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_C2sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_C2sp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_C2sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_C2sp as text
%        str2double(get(hObject,'String')) returns contents of edit_C2sp as a double

sclos2= str2num(get(hObject,'String'));
assignin('base','sclos2',sclos2);
evalin('base','mcr_als.alsOptions.closure.sclos2=sclos2;');
evalin('base','clear sclos2');


% vclos
% oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

% --- Executes on button press in check_Cvar.
function check_Cvar_Callback(hObject, eventdata, handles)
% hObject    handle to check_Cvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_Cvar

vclos=get(hObject,'Value');

if vclos==1;
    iclos=evalin('base','mcr_als.alsOptions.closure.iclos');
    if iclos==1
        set(handles.text_C1,'enable','off');
        set(handles.edit_C1,'enable','off');
        set(handles.text_C1c,'enable','off');
        set(handles.popup_C1c,'enable','off','value',1);
        set(handles.text_C1sp,'enable','on');
        set(handles.edit_C1sp,'enable','on','string',' ');
        set(handles.check_C1sp,'enable','on','value',0);
        set(handles.text_C1v,'enable','on');
        set(handles.edit_C1v,'enable','on','string',' ');
        set(handles.text_C2,'enable','off');
        set(handles.edit_C2,'enable','off');
        set(handles.text_C2c,'enable','off');
        set(handles.popup_C2c,'enable','off','value',1);
        set(handles.text_C2sp,'enable','off');
        set(handles.edit_C2sp,'enable','off','string',' ');
        set(handles.check_C2sp,'enable','off','value',0);
        set(handles.text_C2v,'enable','off');
        set(handles.edit_C2v,'enable','off','string',' ');
        iclos1=1; % has to be the same
        assignin('base','iclos1',iclos1);
        evalin('base','mcr_als.alsOptions.closure.iclos1=iclos1;');
        evalin('base','clear iclos1');
        vc=1; % % has to be the same
        assignin('base','vc',vc);
        evalin('base','mcr_als.alsOptions.closure.vc=vc;');
        evalin('base','clear vc');
        
    elseif iclos==2
        
        
        set(handles.text_C1,'enable','off');
        set(handles.edit_C1,'enable','off');
        set(handles.text_C1c,'enable','off');
        set(handles.popup_C1c,'enable','off','value',1);
        set(handles.text_C1sp,'enable','on');
        set(handles.edit_C1sp,'enable','on','string',' ');
        set(handles.check_C1sp,'enable','on','value',0);
        set(handles.text_C1v,'enable','on');
        set(handles.edit_C1v,'enable','on','string',' ');
        set(handles.text_C2,'enable','off');
        set(handles.edit_C2,'enable','off');
        set(handles.text_C2c,'enable','off');
        set(handles.popup_C2c,'enable','off','value',1);
        set(handles.text_C2sp,'enable','on');
        set(handles.edit_C2sp,'enable','on','string',' ');
        set(handles.check_C2sp,'enable','on','value',0);
        set(handles.text_C2v,'enable','on');
        set(handles.edit_C2v,'enable','on','string',' ');
        iclos1=1; % has to be the same
        assignin('base','iclos1',iclos1);
        evalin('base','mcr_als.alsOptions.closure.iclos1=iclos1;');
        evalin('base','clear iclos1');
        iclos2=1; % has to be the same
        assignin('base','iclos2',iclos2);
        evalin('base','mcr_als.alsOptions.closure.iclos2=iclos2;');
        evalin('base','clear iclos2');
        vc=1; % has to be the same
        assignin('base','vc',vc);
        evalin('base','mcr_als.alsOptions.closure.vc=vc;');
        evalin('base','clear vc');
        
    end
    
elseif vclos==0
    
    iclos=evalin('base','mcr_als.alsOptions.closure.iclos');
    if iclos==1
        set(handles.text_C1,'enable','on');
        set(handles.edit_C1,'enable','on');
        set(handles.text_C1c,'enable','on');
        set(handles.popup_C1c,'enable','on','value',1);
        set(handles.text_C1sp,'enable','on');
        set(handles.edit_C1sp,'enable','on','string',' ');
        set(handles.check_C1sp,'enable','on','value',0);
        set(handles.text_C1v,'enable','off');
        set(handles.edit_C1v,'enable','off','string',' ');
        set(handles.text_C2,'enable','off');
        set(handles.edit_C2,'enable','off');
        set(handles.text_C2c,'enable','off');
        set(handles.popup_C2c,'enable','off','value',1);
        set(handles.text_C2sp,'enable','off');
        set(handles.edit_C2sp,'enable','off','string',' ');
        set(handles.check_C2sp,'enable','off','value',0);
        set(handles.text_C2v,'enable','off');
        set(handles.edit_C2v,'enable','off','string',' ');
        
        iclos1=1; % has to be the same
        assignin('base','iclos1',iclos1);
        evalin('base','mcr_als.alsOptions.closure.iclos1=iclos1;');
        evalin('base','clear iclos1');
        
        vc=1; % has to be the same
        assignin('base','vc',vc);
        evalin('base','mcr_als.alsOptions.closure.vc=vc;');
        evalin('base','clear vc');
        
    elseif iclos==2
                
        set(handles.text_C1,'enable','on');
        set(handles.edit_C1,'enable','on');
        set(handles.text_C1c,'enable','on');
        set(handles.popup_C1c,'enable','on','value',1);
        set(handles.text_C1sp,'enable','on');
        set(handles.edit_C1sp,'enable','on','string',' ');
        set(handles.check_C1sp,'enable','on','value',0);
        set(handles.text_C1v,'enable','off');
        set(handles.edit_C1v,'enable','off','string',' ');
        
        set(handles.text_C2,'enable','on');
        set(handles.edit_C2,'enable','on');
        set(handles.text_C2c,'enable','on');
        set(handles.popup_C2c,'enable','on','value',1);
        set(handles.text_C2sp,'enable','on');
        set(handles.edit_C2sp,'enable','on','string',' ');
        set(handles.check_C2sp,'enable','on','value',0);
        set(handles.text_C2v,'enable','off');
        set(handles.edit_C2v,'enable','off','string',' ');
         
        iclos1=1; % has to be the same 
        assignin('base','iclos1',iclos1);
        evalin('base','mcr_als.alsOptions.closure.iclos1=iclos1;');
        evalin('base','clear iclos1');
        
        iclos2=1; % has to be the same
        assignin('base','iclos2',iclos2);
        evalin('base','mcr_als.alsOptions.closure.iclos2=iclos2;');
        evalin('base','clear iclos2');
        
        vc=1; % has to be the same
        assignin('base','vc',vc);
        evalin('base','mcr_als.alsOptions.closure.vc=vc;');
        evalin('base','clear vc');
        
    end

end;

% --- Executes during object creation, after setting all properties.
function edit_C1v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_C1v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_C1v_Callback(hObject, eventdata, handles)
% hObject    handle to edit_C1v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_C1v as text
%        str2double(get(hObject,'String')) returns contents of edit_C1v as a double

vclos1=get(hObject,'String');
vclos1=evalin('base',vclos1);
assignin('base','vclos1',vclos1);
evalin('base','mcr_als.alsOptions.closure.vclos1=vclos1;');
evalin('base','clear vclos1');


% --- Executes during object creation, after setting all properties.
function edit_C2v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_C2v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_C2v_Callback(hObject, eventdata, handles)
% hObject    handle to edit_C2v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_C2v as text
%        str2double(get(hObject,'String')) returns contents of edit_C2v as a double

vclos2=get(hObject,'String');

vclos2=evalin('base',vclos2);
assignin('base','vclos2',vclos2);
evalin('base','mcr_als.alsOptions.closure.vclos2=vclos2;');
evalin('base','clear vclos2');


% Equality
% *********************************************************************

% --- Executes on button press in check_eqC.
function check_eqC_Callback(hObject, eventdata, handles)
% hObject    handle to check_eqC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_eqC

cselcon=get(hObject,'Value');

if cselcon==1;
    appCorrel=evalin('base','mcr_als.alsOptions.correlation.appCorrelation;');
    if appCorrel==0
        set(handles.text_eqCs,'enable','on');
        set(handles.popup_eqCs,'enable','on','value',1);
        set(handles.text_eqCc,'enable','on');
        set(handles.popup_eqCc,'enable','on','value',1);
    else
        warndlg('no és possible aplicar Equality i Correlation');
        set(hObject,'Value',0);
        cselcon=0;
    end
    
else
    set(handles.text_eqCs,'enable','off');
    set(handles.popup_eqCs,'enable','off','value',1);
    set(handles.text_eqCc,'enable','off');
    set(handles.popup_eqCc,'enable','off','value',1);
         
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''cselcC'');');
    
end;

assignin('base','cselcon',cselcon);
evalin('base','mcr_als.alsOptions.cselcC.cselcon=cselcon;');
evalin('base','clear cselcon');


% --- Executes during object creation, after setting all properties.
function popup_eqCs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_eqCs (see GCBO)
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

% --- Executes on selection change in popup_eqCs.
function popup_eqCs_Callback(hObject, eventdata, handles)
% hObject    handle to popup_eqCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_eqCs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_eqCs

popmenu3=get(handles.popup_eqCs,'String');
pm3=get(handles.popup_eqCs,'Value');

if pm3==1
    
else
    selcsel_char=char([popmenu3(pm3)]);
    selcsel=evalin('base',selcsel_char);
    iisel=find(isfinite(selcsel));
    assignin('base','selcsel',selcsel);
    assignin('base','iisel',iisel);
    assignin('base','selcsel_char',selcsel_char);
    evalin('base','mcr_als.alsOptions.cselcC.iisel=iisel;');
    evalin('base','mcr_als.alsOptions.cselcC.csel=selcsel;');
    evalin('base','mcr_als.alsOptions.cselcC.selcsel_char=selcsel_char;');
    evalin('base','clear selcsel iisel selcsel_char');     
end


% --- Executes during object creation, after setting all properties.
function popup_eqCc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_eqCc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

llista(1)={'select...'};
llista(2)={'equal to'};
llista(3)={'lower than or equal to'};
llista(4)={'greater than'};
set(hObject,'string',llista)

% --- Executes on selection change in popup_eqCc.
function popup_eqCc_Callback(hObject, eventdata, handles)
% hObject    handle to popup_eqCc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_eqCc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_eqCc

type_csel=get(hObject,'Value')-2; 
assignin('base','type_csel',type_csel);
evalin('base','mcr_als.alsOptions.cselcC.type_csel=type_csel;');
evalin('base','clear type_csel');



% Advanced constraints
% *********************************************************************


% Correlation
% ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

% --- Executes on button press in push_correlation.
function push_correlation_Callback(hObject, eventdata, handles)
% hObject    handle to push_correlation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

appCorrel=evalin('base','mcr_als.alsOptions.correlation.appCorrelation;');
appEquality=evalin('base','mcr_als.alsOptions.cselcC.cselcon;');

if appEquality==1
    warndlg('no és possible aplicar Equality i Correlation');
else
    if appCorrel==0
        correlCons;
    else
        correlCons_reload;
    end
end



% Kinetic
% ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

% --- Executes on button press in push_kinetic.
function push_kinetic_Callback(hObject, eventdata, handles)
% hObject    handle to push_kinetic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

kineticHM;



% Continue
% *********************************************************************

% --- Executes on button press in push_continue.
function push_continue_Callback(hObject, eventdata, handles)
% hObject    handle to push_continue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close;
if evalin('base','mcr_als.alsOptions.multi.datamod')==1
    columnModeConstraints;
elseif evalin('base','mcr_als.alsOptions.multi.datamod')==2
    columnMultiConstraints;
elseif evalin('base','mcr_als.alsOptions.multi.datamod')==3
    columnMultiConstraints;
else
       columnModeConstraints;
end

% Reset
% *********************************************************************


% --- Executes on button press in push_reset.
function push_reset_Callback(hObject, eventdata, handles)
% hObject    handle to push_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close;
rowModeConstraints;

% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
