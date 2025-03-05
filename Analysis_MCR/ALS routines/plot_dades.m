function varargout = plot_dades(varargin)
% PLOT_DADES MATLAB code for plot_dades.fig
%      PLOT_DADES, by itself, creates a new PLOT_DADES or raises the existing
%      singleton*.
%
%      H = PLOT_DADES returns the handle to a new PLOT_DADES or the handle to
%      the existing singleton*.
%
%      PLOT_DADES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOT_DADES.M with the given input arguments.
%
%      PLOT_DADES('Property','Value',...) creates a new PLOT_DADES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plot_dades_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plot_dades_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plot_dades

% Last Modified by GUIDE v2.5 09-Feb-2025 22:41:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plot_dades_OpeningFcn, ...
                   'gui_OutputFcn',  @plot_dades_OutputFcn, ...
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


% --- Executes just before plot_dades is made visible.
function plot_dades_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plot_dades (see VARARGIN)

% Choose default command line output for plot_dades
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plot_dades wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plot_dades_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
    
%if matdad_small is in the workspace, only a limited number of spectra will
%be shown 
    if evalin('base','exist(''matdad_small'')')==1
        matdad=evalin('base','matdad_small');
    else
        matdad=evalin('base','mcr_als.data');
    end
    axes(handles.axes1);
    [xi,yi]=size(matdad);
    maxim=max(max(matdad));
    minim=min(min(matdad));
    
    if maxim > 0
        maxim1=maxim+0.2*maxim;
    else
        maxim1=maxim-0.2*abs(maxim)
    end
    
    if minim >= 0
        minim1=minim-0.2*minim;
    else
        minim1=minim-0.2*abs(minim);
    end
    set(gcf,'defaultaxescolororder',jet(xi))
    plot(matdad');axis([1 yi minim1 maxim1]);title('Rows of raw data matrix');xlabel('Channels');ylabel('Intensity (a.u.)');
    zoom on;

% --- Executes during object creation, after setting all properties.
function listbox_sel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox_sel.
function listbox_sel_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_sel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_sel

valor=get(hObject,'Value');

if valor==1
    if evalin('base','exist(''matdad_small'')')==1
        matdad=evalin('base','matdad_small');
    else
        matdad=evalin('base','mcr_als.data');
    end

    axes(handles.axes1);
    [xi,yi]=size(matdad);
    maxim=max(max(matdad));
    minim=min(min(matdad));
    
    if maxim > 0
        maxim1=maxim+0.2*maxim;
    else
        maxim1=maxim-0.2*abs(maxim)
    end
    
    if minim >= 0
        minim1=minim-0.2*minim;
    else
        minim1=minim-0.2*abs(minim);
    end
    set(gcf,'defaultaxescolororder',jet(xi))


    xAxis=get(handles.SetAxis,'Value');
    if xAxis==1
    plot(matdad');axis([1 yi minim1 maxim1]);title('Rows of raw data matrix');xlabel('Channels');ylabel('Intensity (a.u.)');
    zoom on;
    else
    RefAxis=evalin('base','mcr_als.RefAxis');
    plot(RefAxis,matdad');axis([min(RefAxis) max(RefAxis) minim1 maxim1]);title('Rows of raw data matrix');xlabel('Channels');ylabel('Intensity (a.u.)');
    zoom on;
    end

elseif valor==2
    if evalin('base','exist(''matdad_small'')')==1
        matdad=evalin('base','matdad_small');
    else
        matdad=evalin('base','mcr_als.data');
    end
    axes(handles.axes1);
    [xi,yi]=size(matdad);
    maxim=max(max(matdad));
    minim=min(min(matdad));
    
    if maxim > 0
        maxim1=maxim+0.2*maxim;
    else
        maxim1=maxim-0.2*abs(maxim)
    end
    
    if minim >= 0
        minim1=minim-0.2*minim;
    else
        minim1=minim-0.2*abs(minim);
    end
    set(gcf,'defaultaxescolororder',jet(yi))

    xAxis=get(handles.SetAxis,'Value');
    if xAxis==1
    plot(matdad);axis([1 xi minim1 maxim1]);title('Columns of raw data matrix');xlabel('Channels');ylabel('Intensity (a.u.)');
    zoom on;
    else
    RefAxis=evalin('base','mcr_als.RefAxis');
    plot(RefAxis, matdad);axis([min(RefAxis) max(RefAxis) minim1 maxim1]);title('Rows of raw data matrix');xlabel('Channels');ylabel('Intensity (a.u.)');
    zoom on;
    end
end

%NOTE: when the axis is fixed, MATLAB always consider the dataset with the
%proper dimensions. 



% --- Executes on button press in push_close.
function push_close_Callback(hObject, eventdata, handles)
% hObject    handle to push_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close;


% --- Executes on button press in SetAxis.
function SetAxis_Callback(hObject, eventdata, handles)
% hObject    handle to SetAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu

SetAxisMenu=get(handles.SetAxis,'String'); 
xAxis=get(handles.SetAxis,'Value');

if xAxis==1 %the menu label is still selected, so there is no item selection  
    
else
    selec=char([SetAxisMenu(xAxis)]);
    RefAxis=evalin('base',selec);
    assignin('base','RefAxis',RefAxis);
    assignin('base','selec',selec);
    evalin('base','mcr_als.RefAxis=RefAxis;');
    evalin('base','mcr_als.aux.RefAxis_name=selec;');
    % set(handles.text_nom,'string',selec);
    evalin('base','clear RefAxis selec'); 
    
end





% --- Executes during object creation, after setting all properties.
function SetAxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SetAxis (see GCBO)
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


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
