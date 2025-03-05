function varargout = multiwayC(varargin)
% MULTIWAYC MATLAB code for multiwayC.fig
%      MULTIWAYC, by itself, creates a new MULTIWAYC or raises the existing
%      singleton*.
%
%      H = MULTIWAYC returns the handle to a new MULTIWAYC or the handle to
%      the existing singleton*.
%
%      MULTIWAYC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTIWAYC.M with the given input arguments.
%
%      MULTIWAYC('Property','Value',...) creates a new MULTIWAYC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before multiwayC_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to multiwayC_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help multiwayC

% Last Modified by GUIDE v2.5 07-Jan-2014 15:28:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @multiwayC_OpeningFcn, ...
                   'gui_OutputFcn',  @multiwayC_OutputFcn, ...
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


% --- Executes just before multiwayC is made visible.
function multiwayC_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to multiwayC (see VARARGIN)

% Choose default command line output for multiwayC
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes multiwayC wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = multiwayC_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% General Options *********************************************************
% *************************************************************************

% --- Executes during object creation, after setting all properties.
function popup_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

llista(1)={'select type of model...'};
llista(2)={'Non multilinear, only bilinear'}; % iShape = 0
llista(3)={'multilinear, equal shape and synchronization (all species)'}; % iShape = 1
llista(4)={'multilinear without synchronization (all species)'}; % iShape = 2
llista(5)={'multiinear and synchronization for only some species'}; % iShape = 3
llista(6)={'dif. nr. of components in the three modes, interaction models'}; % iShape = 4
set(hObject,'string',llista)

% --- Executes on selection change in popup_type.
function popup_type_Callback(hObject, eventdata, handles)
% hObject    handle to popup_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_type


ishape=get(hObject,'Value')-2;

evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''trilin'');');
multiway_parameters;

assignin('base','ishape',ishape);
evalin('base','mcr_als.alsOptions.trilin.ishape=ishape;');

if ishape==-1 % select option
    set(handles.text_complex,'enable','off');
    set(handles.listbox_complex,'enable','off','value',1);
    set(handles.text_qelem,'enable','off');
    set(handles.listbox_qelem,'enable','off','value',1);
    set(handles.text_qelemNr,'enable','off');
    set(handles.edit_qelemNr,'enable','off','string',' ');
    set(handles.text_Cvector,'enable','off');
    set(handles.edit_Cvector,'enable','off','string',' ');
    set(handles.text_Cinteract,'enable','off');
    set(handles.listbox_Cinteract,'enable','off','value',1);
    set(handles.text_CcompNr,'enable','off');
    set(handles.edit_CcompNr,'enable','off','string',' ');
    set(handles.text_CintP,'enable','off');
    set(handles.edit_CintP,'enable','off','string',' ');
    set(handles.push_done,'enable','off');
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''trilin'');');
    evalin('base','mcr_als.alsOptions.trilin.appTril=0;');
elseif ishape==0
    set(handles.text_complex,'enable','off');
    set(handles.listbox_complex,'enable','off','value',1);
    set(handles.text_qelem,'enable','off');
    set(handles.listbox_qelem,'enable','off','value',1);
    set(handles.text_qelemNr,'enable','off');
    set(handles.edit_qelemNr,'enable','off','string',' ');
    set(handles.text_Cvector,'enable','off');
    set(handles.edit_Cvector,'enable','off','string',' ');
    set(handles.text_Cinteract,'enable','off');
    set(handles.listbox_Cinteract,'enable','off','value',1);
    set(handles.text_CcompNr,'enable','off');
    set(handles.edit_CcompNr,'enable','off','string',' ');
    set(handles.text_CintP,'enable','off');
    set(handles.edit_CintP,'enable','off','string',' ');
    set(handles.push_done,'enable','on');
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''trilin'');');
    % antic
    % evalin('base','mcr_als.alsOptions.trilin.appTril=0;');
    % nou
    evalin('base','mcr_als.alsOptions.trilin.appTril=2;'); % cas bilinear
    evalin('base','mcr_als.alsOptions.trilin.iquadril=0;');
elseif ishape==1
    set(handles.text_complex,'enable','on');
    set(handles.listbox_complex,'enable','on','value',1);
    set(handles.text_qelem,'enable','off');
    set(handles.listbox_qelem,'enable','off','value',1);
    set(handles.text_qelemNr,'enable','off');
    set(handles.edit_qelemNr,'enable','off','string',' ');
    set(handles.text_Cvector,'enable','on');
    
    uns=ones(1,evalin('base','mcr_als.CompNumb.nc;'));
    set(handles.edit_Cvector,'enable','inactive','string',num2str(uns));
    assignin('base','spetric',uns);
    evalin('base','mcr_als.alsOptions.trilin.spetric=spetric;');
    evalin('base','clear spetric');
    
    set(handles.text_Cinteract,'enable','off');
    set(handles.listbox_Cinteract,'enable','off','value',1);
    set(handles.text_CcompNr,'enable','off');
    set(handles.edit_CcompNr,'enable','off','string',' ');
    set(handles.text_CintP,'enable','off');
    set(handles.edit_CintP,'enable','off','string',' ');
    set(handles.push_done,'enable','on');
    % values
    evalin('base','mcr_als.alsOptions.trilin.iquadril=0;');
    evalin('base','mcr_als.alsOptions.trilin.appTril=1;');
elseif ishape==2
    set(handles.text_complex,'enable','on');
    set(handles.listbox_complex,'enable','on','value',1);
    set(handles.text_qelem,'enable','off');
    set(handles.listbox_qelem,'enable','off','value',1);
    set(handles.text_qelemNr,'enable','off');
    set(handles.edit_qelemNr,'enable','off','string',' ');
    set(handles.text_Cvector,'enable','on');

    uns=ones(1,evalin('base','mcr_als.CompNumb.nc;'));
    set(handles.edit_Cvector,'enable','inactive','string',num2str(uns));
    
    assignin('base','spetric',uns);
    evalin('base','mcr_als.alsOptions.trilin.spetric=spetric;');
    evalin('base','clear spetric');
    
    set(handles.text_Cinteract,'enable','off');
    set(handles.listbox_Cinteract,'enable','off','value',1);
    set(handles.text_CcompNr,'enable','off');
    set(handles.edit_CcompNr,'enable','off','string',' ');
    set(handles.text_CintP,'enable','off');
    set(handles.edit_CintP,'enable','off','string',' ');
    set(handles.push_done,'enable','on');
    % values
    evalin('base','mcr_als.alsOptions.trilin.iquadril=0;');
    evalin('base','mcr_als.alsOptions.trilin.appTril=1;');
elseif ishape==3
    set(handles.text_complex,'enable','on');
    set(handles.listbox_complex,'enable','on','value',1);
    set(handles.text_qelem,'enable','off');
    set(handles.listbox_qelem,'enable','off','value',1);
    set(handles.text_qelemNr,'enable','off');
    set(handles.edit_qelemNr,'enable','off','string',' ');
    set(handles.text_Cvector,'enable','on');
    set(handles.edit_Cvector,'enable','on','string',' ');
    set(handles.text_Cinteract,'enable','off');
    set(handles.listbox_Cinteract,'enable','off','value',1);
    set(handles.text_CcompNr,'enable','off');
    set(handles.edit_CcompNr,'enable','off','string',' ');
    set(handles.text_CintP,'enable','off');
    set(handles.edit_CintP,'enable','off','string',' ');
    set(handles.push_done,'enable','off');
    % values
    evalin('base','mcr_als.alsOptions.trilin.iquadril=0;');
    evalin('base','mcr_als.alsOptions.trilin.appTril=1;');
elseif ishape==4
    set(handles.text_complex,'enable','off');
    set(handles.listbox_complex,'enable','off','value',1);
    set(handles.text_qelem,'enable','off');
    set(handles.listbox_qelem,'enable','off','value',1);
    set(handles.text_qelemNr,'enable','off');
    set(handles.edit_qelemNr,'enable','off','string',' ');
    set(handles.text_Cvector,'enable','off');
    set(handles.edit_Cvector,'enable','off','string',' ');
    set(handles.text_Cinteract,'enable','on');
    set(handles.listbox_Cinteract,'enable','on','value',1);
    
    set(handles.text_CcompNr,'enable','on');
    expressionX=['mcr_als.alsOptions.trilin.modeltuckc(1);'];
    valor=evalin('base',expressionX);
    set(handles.edit_CcompNr,'enable','on','string',num2str(valor));
    
    set(handles.text_CintP,'enable','on');
    expressionX=['mcr_als.alsOptions.trilin.spetuckc(1,:);'];
    valor=evalin('base',expressionX);
    set(handles.edit_CintP,'enable','on','string',num2str(valor));
    
    set(handles.push_done,'enable','on');
    % values
    evalin('base','mcr_als.alsOptions.trilin.iquadril=0;');
    evalin('base','mcr_als.alsOptions.trilin.appTril=1;');
end
evalin('base','clear ishape');


% --- Executes during object creation, after setting all properties.
function listbox_complex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_complex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox_complex.
function listbox_complex_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_complex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_complex contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_complex

iquadril=get(hObject,'Value')-1;
assignin('base','iquadril',iquadril);
evalin('base','mcr_als.alsOptions.trilin.iquadril=iquadril;');

if iquadril==1
    set(handles.text_qelem,'enable','on');
    set(handles.listbox_qelem,'enable','on','value',1);
    set(handles.text_qelemNr,'enable','on');
    set(handles.edit_qelemNr,'enable','on','string','0');
elseif iquadril==0
    set(handles.text_qelem,'enable','off');
    set(handles.listbox_qelem,'enable','off','value',1);
    set(handles.text_qelemNr,'enable','off');
    set(handles.edit_qelemNr,'enable','off','string','');    
end
evalin('base','clear iquadril');

% --- Executes during object creation, after setting all properties.
function listbox_qelem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_qelem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox_qelem.
function listbox_qelem_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_qelem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_qelem contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_qelem

Mactiu=get(hObject,'Value');
assignin('base','Mactiu',Mactiu);
evalin('base','mcr_als.alsOptions.trilin.aux.Mactiu=Mactiu;');
% expressionX=['mcr_als.alsOptions.trilin.kinetic.models.model',num2str(modelActiu),'.constants.kactiu=kactiu;'];
% evalin('base',expressionX);
evalin('base','clear Mactiu');

% change the value of the edit box by the value of the memory
expressionX=['mcr_als.alsOptions.trilin.ne(',num2str(Mactiu),');'];
valor=evalin('base',expressionX);
set(handles.edit_qelemNr,'enable','on','string',num2str(valor));


% --- Executes during object creation, after setting all properties.
function edit_qelemNr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_qelemNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_qelemNr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_qelemNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_qelemNr as text
%        str2double(get(hObject,'String')) returns contents of edit_qelemNr as a double

nei=str2num(get(hObject,'String'));
assignin('base','nei',nei);
Mactiu=evalin('base','mcr_als.alsOptions.trilin.aux.Mactiu;');
assignin('base','Mactiu',Mactiu);
expressionX=['mcr_als.alsOptions.trilin.ne(',num2str(Mactiu),')=nei;'];
evalin('base',expressionX);
evalin('base','clear Mactiu nei');


% Concentrations
% ***********************************************************************

% --- Executes during object creation, after setting all properties.
function edit_Cvector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Cvector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Cvector_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Cvector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Cvector as text
%        str2double(get(hObject,'String')) returns contents of edit_Cvector as a double

spetric=str2num(get(hObject,'String'));
assignin('base','spetric',spetric);
evalin('base','mcr_als.alsOptions.trilin.spetric=spetric;');
evalin('base','clear spetric');

set(handles.push_done,'enable','on');



% --- Executes during object creation, after setting all properties.
function listbox_Cinteract_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_Cinteract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_Cinteract.
function listbox_Cinteract_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_Cinteract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_Cinteract contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_Cinteract

Iactiu=get(hObject,'Value');
assignin('base','Iactiu',Iactiu);
evalin('base','mcr_als.alsOptions.trilin.aux.Iactiu=Iactiu;');
% expressionX=['mcr_als.alsOptions.trilin.kinetic.models.model',num2str(modelActiu),'.constants.kactiu=kactiu;'];
% evalin('base',expressionX);
evalin('base','clear Iactiu');

% change the value of the edit box by the value of the memory
expressionX=['mcr_als.alsOptions.trilin.modeltuckc(',num2str(Iactiu),');'];
valor=evalin('base',expressionX);
set(handles.edit_CcompNr,'enable','on','string',num2str(valor));

expressionX=['mcr_als.alsOptions.trilin.spetuckc(',num2str(Iactiu),',:);'];
valor=evalin('base',expressionX);
set(handles.edit_CintP,'enable','on','string',num2str(valor));


% --- Executes during object creation, after setting all properties.
function edit_CcompNr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CcompNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_CcompNr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CcompNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CcompNr as text
%        str2double(get(hObject,'String')) returns contents of edit_CcompNr as a double

Ccomp=str2num(get(hObject,'String'));
assignin('base','Ccomp',Ccomp);
Iactiu=evalin('base','mcr_als.alsOptions.trilin.aux.Iactiu;');
assignin('base','Iactiu',Iactiu);
expressionX=['mcr_als.alsOptions.trilin.modeltuckc(',num2str(Iactiu),')=Ccomp;'];
evalin('base',expressionX);
evalin('base','clear Iactiu Ccomp');


% --- Executes during object creation, after setting all properties.
function edit_CintP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CintP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_CintP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CintP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CintP as text
%        str2double(get(hObject,'String')) returns contents of edit_CintP as a double

Ccomp=str2num(get(hObject,'String'));
assignin('base','Ccomp',Ccomp);
Iactiu=evalin('base','mcr_als.alsOptions.trilin.aux.Iactiu;');
assignin('base','Iactiu',Iactiu);
expressionX=['mcr_als.alsOptions.trilin.spetuckc(',num2str(Iactiu),',:)=Ccomp;'];
evalin('base',expressionX);
evalin('base','clear Iactiu Ccomp');




% Spectra
% ***********************************************************************

% not available


% buttons
% ************************************************************


% --- Executes on button press in push_back.
function push_back_Callback(hObject, eventdata, handles)
% hObject    handle to push_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''trilin'');');
    evalin('base','mcr_als.alsOptions.trilin.appTril=0;');
    close;

% --- Executes on button press in push_done.
function push_done_Callback(hObject, eventdata, handles)
% hObject    handle to push_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;
