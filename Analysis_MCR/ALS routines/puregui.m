function varargout = puregui(varargin)
% PUREGUI M-file for puregui.fig
%      PUREGUI, by itself, creates a new PUREGUI or raises the existing
%      singleton*.
%
%      H = PUREGUI returns the handle to a new PUREGUI or the handle to
%      the existing singleton*.
%
%      PUREGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PUREGUI.M with the given input arguments.
%
%      PUREGUI('Property','Value',...) creates a new PUREGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before puregui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to puregui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help puregui

% Last Modified by GUIDE v2.5 03-Oct-2010 23:30:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @puregui_OpeningFcn, ...
                   'gui_OutputFcn',  @puregui_OutputFcn, ...
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


% --- Executes just before puregui is made visible.
function puregui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to puregui (see VARARGIN)

% Choose default command line output for puregui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes puregui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = puregui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% initialization
evalin('base','mcr_als.InitEstim.Sort=0;');


% ***********************************************************************
% Edit noise

% --- Executes during object creation, after setting all properties.
function edit_noise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

evalin('base','mcr_als.InitEstim.pureNoiseLevel=10;');

function edit_noise_Callback(hObject, eventdata, handles)
% hObject    handle to edit_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_noise as text
%        str2double(get(hObject,'String')) returns contents of edit_noise as a double

noiseLevel=str2double(get(hObject,'String'));
assignin('base','noiseLevel',noiseLevel);
evalin('base','mcr_als.InitEstim.pureNoiseLevel=noiseLevel;');
evalin('base','clear noiseLevel');


% ***********************************************************************
% List of variables

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% ***********************************************************************
% PopUp directions

% --- Executes during object creation, after setting all properties.
function popup_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

llista(1)={'select...'};
llista(2)={'Concentrations'};
llista(3)={'Spectra'};
set(hObject,'string',llista)

% --- Executes on selection change in popup_dir.
function popup_dir_Callback(hObject, eventdata, handles)
% hObject    handle to popup_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_dir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_dir

direction=get(hObject,'Value')-1; 

if direction==0
    set(handles.push_fer,'enable','off');
    set(handles.push_OK,'enable','off');
else
    set(handles.push_fer,'enable','on');
end

assignin('base','direction',direction);
evalin('base','mcr_als.InitEstim.pureDirection=direction;');
evalin('base','clear direction');


% ***********************************************************************
% Exe push buttons

% --- Executes on button press in push_fer.
function push_fer_Callback(hObject, eventdata, handles)
% hObject    handle to push_fer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d=evalin('base','mcr_als.data;');
nr=evalin('base','mcr_als.CompNumb.nc;');
if (evalin('base','mcr_als.InitEstim.pureDirection;')==2)
    d=d';
end
f=evalin('base','mcr_als.InitEstim.pureNoiseLevel;');
[nrow,ncol]=size(d);

% calculation of the purity spectrum

f=f/100;
s=std(d);
m=mean(d);
ll=s.*s+m.*m;
f=max(m)*f;
p=s./(m+f);

[mp,imp(1)]=max(p);

l=sqrt((s.*s+(m+f).*(m+f)));

for j=1:ncol,
    dl(:,j)=d(:,j)./l(j);
end
c=(dl'*dl)./nrow;

% calculation of the weights
% first weight

w(1,:)=ll./(l.*l);
p(1,:)=w(1,:).*p(1,:);
s(1,:)=w(1,:).*s(1,:);

% next weights

for i=2:nr,
    for j=1:ncol,
        [dm]=wmat(c,imp,i,j);
        w(i,j)=det(dm);
        p(i,j)=p(1,j).*w(i,j);
        s(i,j)=s(1,j).*w(i,j);
    end

    % next purest and standard deviation spectrum

    [mp(i),imp(i)]=max(p(i,:));
%     disp('next purest variable: ');disp(imp(i))
end

for i=1:nr,
    impi=imp(i);
    sp(1:nrow,i)=d(1:nrow,impi);
end

sp=normv2(sp');
axes(handles.axes1)
plot(sp');

for i=1:length(imp)
    Scell(i)={imp(i)};
end

% list of purest variables
assignin('base','imp',imp);
evalin('base','mcr_als.InitEstim.indices=imp;');
evalin('base','clear imp');

%set(handles.list_eigval,'string',Scell)
set(handles.listbox1,'string',Scell,'value',1);
set(handles.push_OK,'enable','on');

assignin('base','sp',sp);
evalin('base','mcr_als.InitEstim.iniesta=sp;');
evalin('base','clear sp');


% ***********************************************************************
% Options buttons

% --- Executes on button press in push_OK.
function push_OK_Callback(hObject, eventdata, handles)
% hObject    handle to push_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evalin('base','mcr_als.aux.estat=2;');
evalin('base','mcr_als.InitEstim.method=''Pur'';');

% addition for sorting
if evalin('base','mcr_als.InitEstim.Sort') == 0;
    % nothing to do
else
    % sort
    indexs=evalin('base','mcr_als.InitEstim.indices;');
    sp=evalin('base','mcr_als.InitEstim.iniesta;');
    [x,i]=sort(indexs);
    sp=sp(i,:);   
    
    assignin('base','sp',sp);
    evalin('base','mcr_als.InitEstim.iniesta=sp;');
    evalin('base','clear sp');
    
end

close;


% --- Executes on button press in push_Print.
function push_Print_Callback(hObject, eventdata, handles)
% hObject    handle to push_Print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(gcf,'PaperPositionMode','auto');
print('PUREGUI','-dpdf');

% --- Executes on button press in push_Copy.
function push_Copy_Callback(hObject, eventdata, handles)
% hObject    handle to push_Copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

print -dbitmap;

% --- Executes on button press in push_canc.
function push_canc_Callback(hObject, eventdata, handles)
% hObject    handle to push_canc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evalin('base','mcr_als.aux.estat=3;');
evalin('base','mcr_als=rmfield(mcr_als,''InitEstim'');');
close;


% --- Executes on button press in check_sort.
function check_sort_Callback(hObject, eventdata, handles)
% hObject    handle to check_sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_sort

pon=get(handles.check_sort,'Value');

if pon==1; % yes
evalin('base','mcr_als.InitEstim.Sort=1;');
else % no
evalin('base','mcr_als.InitEstim.Sort=0;');
end
