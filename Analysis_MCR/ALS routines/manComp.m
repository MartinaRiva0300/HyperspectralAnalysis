function varargout = manComp(varargin)
% MANCOMP M-file for manComp.fig
%      MANCOMP, by itself, creates a new MANCOMP or raises the existing
%      singleton*.
%
%      H = MANCOMP returns the handle to a new MANCOMP or the handle to
%      the existing singleton*.
%
%      MANCOMP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANCOMP.M with the given input arguments.
%
%      MANCOMP('Property','Value',...) creates a new MANCOMP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before manComp_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to manComp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help manComp

% Last Modified by GUIDE v2.5 29-Jan-2008 13:43:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @manComp_OpeningFcn, ...
                   'gui_OutputFcn',  @manComp_OutputFcn, ...
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


% --- Executes just before manComp is made visible.
function manComp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manComp (see VARARGIN)

% Choose default command line output for manComp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes manComp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = manComp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ***********************************************************
% ***********************************************************
% ***********************************************************

% --- Executes during object creation, after setting all properties.
function edit_NC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_NC_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NC as text
%        str2double(get(hObject,'String')) returns contents of edit_NC as a double

% @post. the user writes the NC -> continue -> mcr_als and election

nc=str2num(get(hObject,'String'));

if nc>0
assignin('base','nc',nc);
evalin('base','mcr_als.CompNumb.nc=nc;');
evalin('base','mcr_als.CompNumb.method=''manual'';');
evalin('base','clear nc');
set(handles.push_Cont,'enable','on');
else
set(handles.push_Cont,'enable','off'); 
if ((evalin('base','exist(''mcr_als.CompNumb'')'))==1);
    evalin('base','mcr_als=rmfield(mcr_als,''CompNumb'');');end;
end

% ***********************************************************
% ***********************************************************
% ***********************************************************

% --- Executes on button press in push_Cont.
function push_Cont_Callback(hObject, eventdata, handles)
% hObject    handle to push_Cont (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evalin('base','mcr_als.aux.estat=2;');
close;


% ***********************************************************
% ***********************************************************
% ***********************************************************

% --- Executes on button press in push_Quit.
function push_Quit_Callback(hObject, eventdata, handles)
% hObject    handle to push_Quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evalin('base','mcr_als.aux.estat=3;');
evalin('base','mcr_als=rmfield(mcr_als,''CompNumb'');');
close;
% instruction to delete the value of the determination of NC is required

