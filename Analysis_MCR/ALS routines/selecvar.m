function varargout = selecvar(varargin)
% SELECVAR M-file for selecvar.fig
%      SELECVAR, by itself, creates a new SELECVAR or raises the existing
%      singleton*.
%
%      H = SELECVAR returns the handle to a new SELECVAR or the handle to
%      the existing singleton*.
%
%      SELECVAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECVAR.M with the given input arguments.
%
%      SELECVAR('Property','Value',...) creates a new SELECVAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before selecvar_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to selecvar_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help selecvar

% Last Modified by GUIDE v2.5 07-Feb-2008 17:50:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @selecvar_OpeningFcn, ...
                   'gui_OutputFcn',  @selecvar_OutputFcn, ...
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


% --- Executes just before selecvar is made visible.
function selecvar_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to selecvar (see VARARGIN)

% Choose default command line output for selecvar
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes selecvar wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = selecvar_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

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

Sd=evalin('base','mcr_als.aux.a');

for i=1:length(Sd)
    Scell(i)={Sd(i)};
end

%set(handles.list_eigval,'string',Scell)
set(hObject,'string',Sd,'value',1)

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


Number=get(hObject,'Value');

assignin('base','Number',Number);
evalin('base','mcr_als.aux.Number=Number;');
evalin('base','clear Number');

% --- Executes on button press in push_OK.
function push_OK_Callback(hObject, eventdata, handles)
% hObject    handle to push_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evalin('base','mcr_als.aux.estat=0;');
close;

% --- Executes on button press in push_cancel.
function push_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to push_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evalin('base','mcr_als.aux.estat=0;');
close;
