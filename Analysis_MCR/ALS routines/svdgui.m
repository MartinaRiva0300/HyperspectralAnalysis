function varargout = svdgui(varargin)
% SVDGUI M-file for svdgui.fig
%      SVDGUI, by itself, creates a new SVDGUI or raises the existing
%      singleton*.
%
%      H = SVDGUI returns the handle to a new SVDGUI or the handle to
%      the existing singleton*.
%
%      SVDGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SVDGUI.M with the given input arguments.
%
%      SVDGUI('Property','Value',...) creates a new SVDGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before svdgui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to svdgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help svdgui

% Last Modified by GUIDE v2.5 04-Feb-2008 12:58:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @svdgui_OpeningFcn, ...
    'gui_OutputFcn',  @svdgui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before svdgui is made visible.
function svdgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to svdgui (see VARARGIN)



% Choose default command line output for svdgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% matriu_svd=evalin('base','mcr_als.data');
% [a,b]=size(matriu_svd);
% evalin('base','mcr_als.aux.estat_svd=1;');
% status=1;
% 
% if (a>10000 | b>10000)
%     svd_selection;
%     while evalin('base','mcr_als.aux.estat_svd')==1
% pause(0.3);
%     end
%     nrsvd=evalin('base','mcr_als.aux.nrsvd;');
%     [U,S,V]=svds(matriu_svd,nrsvd);
% else
%     [U,S,V]=svd(matriu_svd,'econ');
% end
axes(handles.svalues);
S=evalin('base','mcr_als.CompNumb.S;');
plot(diag(S),'o');
zoom on;

% UIWAIT makes svdgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = svdgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Eigenvalues list
% *************************************************************************

% --- Executes during object creation, after setting all properties.
function list_eigval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_eigval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% matriu_svd=evalin('base','mcr_als.data');
% [U,S,V]=svd(matriu_svd);
% assignin('base','U',U);
% assignin('base','S',S);
% assignin('base','V',V);
% evalin('base','mcr_als.CompNumb.U=U;');
% evalin('base','mcr_als.CompNumb.S=S;');
% evalin('base','mcr_als.CompNumb.V=V;');
% evalin('base','clear U S V');

S=evalin('base','mcr_als.CompNumb.S');
Sd=diag(S);
for i=1:length(Sd)
    Scell(i)={Sd(i)};
end

set(hObject,'string',Scell,'value',1)


% --- Executes on selection change in list_eigval.
function list_eigval_Callback(hObject, eventdata, handles)
% hObject    handle to list_eigval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_eigval contents as cell array
%        contents{} returns selected item from list_eigval

U=evalin('base','mcr_als.CompNumb.U');
S=evalin('base','mcr_als.CompNumb.S');
ns=get(hObject,'Value');
axes(handles.svectors);
plot(U(:,1:ns));
zoom on;

normx=get(handles.check_norm,'Value');
x=diag(S);
axes(handles.svalues);
if normx==1;
plot(x/max(x),'o');
elseif normx==0;
    plot(x,'o');
end
zoom on;

set(handles.edit_nc,'string',num2str(ns));
set(handles.push_OK,'enable','on');
assignin('base','nc',ns);
evalin('base','mcr_als.CompNumb.nc=nc;');
evalin('base','clear nc');

% --- Executes on button press in check_norm.
function check_norm_Callback(hObject, eventdata, handles)
% hObject    handle to check_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_norm

pon=get(handles.check_norm,'Value');

if pon==1;
    axes(handles.svalues);
    evalin('base',' plot(diag(mcr_als.CompNumb.S)/max(diag(mcr_als.CompNumb.S)),''o'')');
    zoom on;
else
    axes(handles.svalues);
    evalin('base',' plot(diag(mcr_als.CompNumb.S),''o'')');
    zoom on;
end

% Opcions de bigsvd
% ************************************************************************
% --- Executes on button press in check_bsvd.
%function check_bsvd_Callback(hObject, eventdata, handles)
% hObject    handle to check_bsvd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_bsvd

%bigpca=get(handles.check_bsvd,'Value');


% Butons
% ************************************************************************
% --- Executes on button press in push_OK.
function push_OK_Callback(hObject, eventdata, handles)
% hObject    handle to push_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evalin('base','mcr_als.aux.estat=2;');
evalin('base','mcr_als.CompNumb.method=''SVD'';');
close;

% --- Executes on button press in push_canc.
function push_canc_Callback(hObject, eventdata, handles)
% hObject    handle to push_canc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evalin('base','mcr_als.aux.estat=3;');
evalin('base','mcr_als=rmfield(mcr_als,''CompNumb'');');
close;

% --- Executes on button press in push_copy.
function push_copy_Callback(hObject, eventdata, handles)
% hObject    handle to push_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

print -dbitmap;


% --- Executes on button press in push_print.
function push_print_Callback(hObject, eventdata, handles)
% hObject    handle to push_print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(gcf,'PaperPositionMode','auto');
print('SVDGUI','-dpdf');


% --- Executes when figure1 window is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edit_nc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nc as text
%        str2double(get(hObject,'String')) returns contents of edit_nc as a double

% --- Executes during object creation, after setting all properties.
function edit_nc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




