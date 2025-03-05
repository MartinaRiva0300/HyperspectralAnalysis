function varargout = manEst(varargin)
% MANEST M-file for manEst.fig
%      MANEST, by itself, creates a new MANEST or raises the existing
%      singleton*.
%
%      H = MANEST returns the handle to a new MANEST or the handle to
%      the existing singleton*.
%
%      MANEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANEST.M with the given input arguments.
%
%      MANEST('Property','Value',...) creates a new MANEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before manEst_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to manEst_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help manEst

% Last Modified by GUIDE v2.5 30-Jan-2008 14:21:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @manEst_OpeningFcn, ...
                   'gui_OutputFcn',  @manEst_OutputFcn, ...
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


% --- Executes just before manEst is made visible.
function manEst_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manEst (see VARARGIN)

% Choose default command line output for manEst
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes manEst wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = manEst_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% variable selection
% *************************************************************************

% --- Executes during object creation, after setting all properties.
function pop_llista_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_llista (see GCBO)
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
lsv(1)={'select a variable...'};
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


% --- Executes on selection change in pop_llista.
function pop_llista_Callback(hObject, eventdata, handles)
% hObject    handle to pop_llista (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_llista contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_llista

popmenu2=get(handles.pop_llista,'String');
pm2=get(handles.pop_llista,'Value');

if pm2==1

else
    selec2=char([popmenu2(pm2)]);
    iniesta=evalin('base',selec2);

    axes(handles.axes_plot);
    [xi,yi]=size(iniesta);
    
    maxim=max(max(iniesta));
    minim=min(min(iniesta));
    
    if maxim > 0
        maxim1=maxim+0.2*maxim;
    else
        maxim1=maxim-0.2*abs(maxim)
    end
    
    if minim >= 0
        minim1=minim+0.2*minim;
    else
        minim1=minim-0.2*abs(minim);
    end

    data=evalin('base','mcr_als.data;');
    [md,nd]=size(data);

    if nd==yi
        plot(iniesta');
        axis([1 yi minim1 maxim1]);
    elseif md==xi
        plot(iniesta);
        axis([1 xi minim1 maxim1]);
    end
        
    
    zoom on;

    assignin('base','iniesta',iniesta);
    evalin('base','mcr_als.InitEstim.iniesta=iniesta;');
    evalin('base','clear iniesta');  
        
end


% continue
% *************************************************************************

% --- Executes on button press in but_confi.
function but_confi_Callback(hObject, eventdata, handles)
% hObject    handle to but_confi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

NC=evalin('base','mcr_als.CompNumb.nc;');
NCest=min(size(evalin('base','mcr_als.InitEstim.iniesta')));

if NC==NCest
    evalin('base','mcr_als.aux.estat=2;');
    evalin('base','mcr_als.InitEstim.method=''man'';');
    close;
else
    warndlg('The number of components of the initial estimation and the previously selected does not match');
end
    

% exit
% *************************************************************************

% --- Executes on button press in but_surt.
function but_surt_Callback(hObject, eventdata, handles)
% hObject    handle to but_surt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evalin('base','mcr_als.aux.estat=3;');
evalin('base','mcr_als=rmfield(mcr_als,''InitEstim'');');
close;


