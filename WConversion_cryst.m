function varargout = WConversion_cryst(varargin)
% WCONVERSION_CRYST MATLAB code for WConversion_cryst.fig
%      WCONVERSION_CRYST, by itself, creates a new WCONVERSION_CRYST or raises the existing
%      singleton*.
%
%      H = WCONVERSION_CRYST returns the handle to a new WCONVERSION_CRYST or the handle to
%      the existing singleton*.
%
%      WCONVERSION_CRYST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WCONVERSION_CRYST.M with the given input arguments.
%
%      WCONVERSION_CRYST('Property','Value',...) creates a new WCONVERSION_CRYST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WConversion_cryst_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WConversion_cryst_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WConversion_cryst

% Last Modified by GUIDE v2.5 23-Nov-2022 19:34:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WConversion_cryst_OpeningFcn, ...
                   'gui_OutputFcn',  @WConversion_cryst_OutputFcn, ...
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


% --- Executes just before WConversion_cryst is made visible.
function WConversion_cryst_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WConversion_cryst (see VARARGIN)

global c laser laser_frequency wavelength frequency wavenumber  ...
       pseudofrequency pseudowavelength step delay nyquist ...
       spectralResolution apexAngle temperature crystal Delta_n ...
       calibration_button f_0 fr_real0 calibration_derivative ...
       apodization_type Delta_f

% Choose default command line output for WConversion_cryst
handles.output = hObject;

% Update handles structure
c=299792458; %speed of light [m/s]

%initialization of the pop-up menu
set(handles.crystalMenu,'Value',1);
set(handles.nyquistMenu,'Value',1);
set(handles.spResolutionMenu,'Value',1);
set(handles.apodizationMenu,'Value',1);

%initialization of the calibration button
set(handles.calibration_enable,'Value',0);

calibration_button=0; %no calibration loaded
set(handles.calibration, 'String', '(none)'); %text in the GUI that indicates loaded calibration
f_0=[]; fr_real0=[]; calibration_derivative=[]; %no calibration axis loaded

temperature=23; %temperature initialization [°C]
apexAngle=10; %wedge apex angle initialization [deg.]
laser=532; %laser wavelength initialized [nm]
laser_frequency=c/laser*1e-3;
wavelength=600; %wavelength initialized [nm]
frequency=c/wavelength*1e-3; %frequency initialization calculated from wavelength [THz]
wavenumber=(laser_frequency-frequency)/c*1e10; %wavenumber initialization calculated from other parameters [cm-1]

switch get(handles.crystalMenu,'Value')
    case 1 %YVO4
        crystal=23; %number for YVO4 in nCryst function
    case 2 %alpha-BBO
        crystal=42; %number for alpha-BBO in nCryst function
    case 3 %Calomel
        crystal=47; %number of Calomel in nCryst function
end

Delta_n=abs(diff(nCryst(crystal,wavelength.*1e-3,temperature)));

pseudowavelength=(1/frequency)/(Delta_n/c*sin(apexAngle/180*pi))*1e-6; %pseudowavelength initialization (optical cycle in term of wedge translation) [um]
pseudofrequency=pseudowavelength^(-1); %pseudofrequency initialization [um-1]

delay=5; %delay initialization [fs]
step=delay/(Delta_n/c*sin(apexAngle/180*pi))*1e-9; %step initialization calculated from the delay [um]

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

switch get(handles.apodizationMenu,'Value')
    case 1
        apodization_type=0; %rectangular
    case 2
        apodization_type=1; %Happ-Genzel
    case 3
        apodization_type=2; %3-term Blackmann-Harris
    case 4
        apodization_type=3; %4-term Blackmann-Harris
    case 5
        apodization_type=4; %Triangular
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

set(handles.apexAngleEdit,'String',num2str(apexAngle));
set(handles.temperatureEdit,'String',num2str(temperature));
set(handles.wavelengthEdit,'String',num2str(wavelength));
set(handles.frequencyEdit,'String',num2str(frequency));
set(handles.wavenumberEdit,'String',num2str(wavenumber));
set(handles.laserEdit,'String',num2str(laser));
set(handles.pseudowavelengthEdit,'String',num2str(pseudowavelength));
set(handles.pseudofrequencyEdit,'String',num2str(pseudofrequency));
set(handles.delayEdit,'String',num2str(delay));
set(handles.stepEdit,'String',num2str(step));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


guidata(hObject, handles);

% UIWAIT makes WConversion_cryst wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WConversion_cryst_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in crystalMenu.
function crystalMenu_Callback(hObject, eventdata, handles)
% hObject    handle to crystalMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature nyquist spectralResolution ...
       apodization_type
% Hints: contents = cellstr(get(hObject,'String')) returns crystalMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from crystalMenu
switch get(hObject,'Value')
    case 1 %YVO4
        crystal=23; %number of YVO4 in nCryst function
    case 2 %alpha-BBO
        crystal=42; %number of alpha-BBO in nCryst function
    case 3 %Calomel
        crystal=47; %number of Calomel in nCryst function
end

Delta_n=abs(diff(nCryst(crystal,wavelength.*1e-3,temperature)));

pseudowavelength=(1/frequency)/(Delta_n/c*sin(apexAngle/180*pi))*1e-6; %optical cycle [um]
pseudofrequency=pseudowavelength^(-1); %[um-1]

delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.pseudowavelengthEdit,'String',num2str(pseudowavelength));
set(handles.pseudofrequencyEdit,'String',num2str(pseudofrequency));
set(handles.delayEdit,'String',num2str(delay));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));




% --- Executes during object creation, after setting all properties.
function crystalMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crystalMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function apexAngleEdit_Callback(hObject, eventdata, handles)
% hObject    handle to apexAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle nyquist spectralResolution ...
       apodization_type
% Hints: get(hObject,'String') returns contents of apexAngleEdit as text
%        str2double(get(hObject,'String')) returns contents of apexAngleEdit as a double

%update the Edit window
apexAngle=str2num(get(handles.apexAngleEdit,'String'));
set(handles.apexAngleEdit,'String',apexAngle);

pseudowavelength=(1/frequency)/(Delta_n/c*sin(apexAngle/180*pi))*1e-6; %optical cycle [um]
pseudofrequency=pseudowavelength^(-1); %[um-1]

delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.pseudowavelengthEdit,'String',num2str(pseudowavelength));
set(handles.pseudofrequencyEdit,'String',num2str(pseudofrequency));
set(handles.delayEdit,'String',num2str(delay));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


% --- Executes during object creation, after setting all properties.
function apexAngleEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apexAngleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wavelengthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to wavelengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature laser_frequency ...
       wavenumber nyquist spectralResolution ...
       calibration_button f_0 fr_real0 calibration_derivative ...
       apodization_type
% Hints: get(hObject,'String') returns contents of wavelengthEdit as text
%        str2double(get(hObject,'String')) returns contents of wavelengthEdit as a double
   
%update the Edit window
wavelength=str2num(get(handles.wavelengthEdit,'String'));
set(handles.wavelengthEdit,'String',wavelength);

frequency=c/wavelength*1e-3; %[THz]
wavenumber=(laser_frequency-frequency)/c*1e10; %[cm-1]

if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
    Delta_n=abs(diff(nCryst(crystal,wavelength.*1e-3,temperature)));
    pseudowavelength=(1/frequency)/(Delta_n/c*sin(apexAngle/180*pi))*1e-6; %optical cycle [um]
    pseudofrequency=pseudowavelength^(-1); %[um-1]
    
    delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]
    
else %if calibration is loaded: calculation using experimental TWINS calibration
    pseudofrequency=interp1(fr_real0,f_0,frequency); %[um-1]
    pseudowavelength=pseudofrequency^(-1); %[um-1]
    
    delay=interp1(f_0,calibration_derivative,pseudofrequency).^(-1).*step.*1e3; %[fs] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
end

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.frequencyEdit,'String',num2str(frequency));
set(handles.wavenumberEdit,'String',num2str(wavenumber));
set(handles.pseudowavelengthEdit,'String',num2str(pseudowavelength));
set(handles.pseudofrequencyEdit,'String',num2str(pseudofrequency));
set(handles.delayEdit,'String',num2str(delay));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));

% --- Executes during object creation, after setting all properties.
function wavelengthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wavelengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frequencyEdit_Callback(hObject, eventdata, handles)
% hObject    handle to frequencyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature laser_frequency ...
       wavenumber nyquist spectralResolution ...
       calibration_button fr_real0 f_0 calibration_derivative ...
       apodization_type
% Hints: get(hObject,'String') returns contents of frequencyEdit as text
%        str2double(get(hObject,'String')) returns contents of frequencyEdit as a double
   
%update the Edit window
frequency=str2num(get(handles.frequencyEdit,'String'));
set(handles.frequencyEdit,'String',frequency);

wavelength=c/frequency*1e-3; %[nm]
wavenumber=(laser_frequency-frequency)/c*1e10; %[cm-1]

if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
    Delta_n=abs(diff(nCryst(crystal,wavelength.*1e-3,temperature)));
    pseudowavelength=(1/frequency)/(Delta_n/c*sin(apexAngle/180*pi))*1e-6; %optical cycle [um]
    pseudofrequency=pseudowavelength^(-1); %[um-1]
    
    delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]
    
else %if calibration is loaded: calculation using experimental TWINS calibration
    pseudofrequency=interp1(fr_real0,f_0,frequency); %[um-1]
    pseudowavelength=pseudofrequency^(-1); %[um-1]
    
    delay=interp1(f_0,calibration_derivative,pseudofrequency).^(-1).*step.*1e3; %[fs] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
end

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.wavelengthEdit,'String',num2str(wavelength));
set(handles.wavenumberEdit,'String',num2str(wavenumber));
set(handles.pseudowavelengthEdit,'String',num2str(pseudowavelength));
set(handles.pseudofrequencyEdit,'String',num2str(pseudofrequency));
set(handles.delayEdit,'String',num2str(delay));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


% --- Executes during object creation, after setting all properties.
function frequencyEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frequencyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wavenumberEdit_Callback(hObject, eventdata, handles)
% hObject    handle to wavenumberEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature laser_frequency ...
       wavenumber nyquist spectralResolution ...
       calibration_button fr_real0 f_0 calibration_derivative ...
       apodization_type
% Hints: get(hObject,'String') returns contents of wavenumberEdit as text
%        str2double(get(hObject,'String')) returns contents of wavenumberEdit as a double
   
%update the Edit window
wavenumber=str2num(get(handles.wavenumberEdit,'String'));
set(handles.wavenumberEdit,'String',wavenumber);

frequency=laser_frequency-wavenumber*c*1e-10; %[THz]
wavelength=c/frequency*1e-3; %[nm]

if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
    Delta_n=abs(diff(nCryst(crystal,wavelength.*1e-3,temperature)));
    pseudowavelength=(1/frequency)/(Delta_n/c*sin(apexAngle/180*pi))*1e-6; %optical cycle [um]
    pseudofrequency=pseudowavelength^(-1); %[um-1]
    
    delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]
    
else %if calibration is loaded: calculation using experimental TWINS calibration
    pseudofrequency=interp1(fr_real0,f_0,frequency); %[um-1]
    pseudowavelength=pseudofrequency^(-1); %[um-1]
    
    delay=interp1(f_0,calibration_derivative,pseudofrequency).^(-1).*step.*1e3; %[fs] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
end

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.wavelengthEdit,'String',num2str(wavelength));
set(handles.frequencyEdit,'String',num2str(frequency));
set(handles.pseudowavelengthEdit,'String',num2str(pseudowavelength));
set(handles.pseudofrequencyEdit,'String',num2str(pseudofrequency));
set(handles.delayEdit,'String',num2str(delay));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


% --- Executes during object creation, after setting all properties.
function wavenumberEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wavenumberEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function laserEdit_Callback(hObject, eventdata, handles)
% hObject    handle to laserEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature laser_frequency ...
       wavenumber laser
% Hints: get(hObject,'String') returns contents of laserEdit as text
%        str2double(get(hObject,'String')) returns contents of laserEdit as a double
   
%update the Edit window
laser=str2num(get(handles.laserEdit,'String'));
set(handles.laserEdit,'String',laser);

laser_frequency=c/laser*1e-3; %[THz]
wavenumber=(laser_frequency-frequency)/c*1e10; %[cm-1]

%update changed parameters
set(handles.wavenumberEdit,'String',num2str(wavenumber));


% --- Executes during object creation, after setting all properties.
function laserEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to laserEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pseudowavelengthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pseudowavelengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature laser_frequency ...
       wavenumber nyquist spectralResolution ...
       calibration_button fr_real0 f_0 calibration_derivative ...
       apodization_type
% Hints: get(hObject,'String') returns contents of pseudowavelengthEdit as text
%        str2double(get(hObject,'String')) returns contents of pseudowavelengthEdit as a double
   
%update the Edit window
pseudowavelength=str2num(get(handles.pseudowavelengthEdit,'String'));
set(handles.pseudowavelengthEdit,'String',pseudowavelength);

pseudofrequency=pseudowavelength^(-1); %[um-1]

if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
    %WAVELENGTH AS FUNCTION OF THE PSEUDOWAVELENGTH (OPTICAL CYCLE IN um)
    %wavelength=Delta_n*sin(apexAngle)*pseudowavelength
    %as Delta_n depends on wavelength a numerical method to solve this equation
    %in wavelength is necessary
    %We solve f(wavelength)=0 with the Newton method
    f=@(variable) variable-abs(diff(nCryst(crystal,variable*1e-3,temperature))).*sin(apexAngle/180*pi).*pseudowavelength.*1e3; %[nm]
    initialGuess=Delta_n.*sin(apexAngle/180*pi).*pseudowavelength.*1e3; %initial guess [nm] exploiting the actual Delta_n
    iterations=100; %maximum number of iterations in the method
    tollerance=1e-5; %error tollerance [nm]
    wavelength=newton_method(f,initialGuess,iterations,tollerance); %[nm]

    frequency=c/wavelength*1e-3; %[THz]
    
    Delta_n=abs(diff(nCryst(crystal,wavelength.*1e-3,temperature)));

    delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]
    
else %if calibration is loaded: calculation using experimental TWINS calibration
    frequency=interp1(f_0,fr_real0,pseudofrequency); %[THz]
    wavelength=c/frequency*1e-3; %[nm]
    
    delay=interp1(f_0,calibration_derivative,pseudofrequency).^(-1).*step.*1e3; %[fs] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
end


wavenumber=(laser_frequency-frequency)/c*1e10; %[cm-1]


switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.wavelengthEdit,'String',num2str(wavelength));
set(handles.frequencyEdit,'String',num2str(frequency));
set(handles.wavenumberEdit,'String',num2str(wavenumber));
set(handles.pseudofrequencyEdit,'String',num2str(pseudofrequency));
set(handles.delayEdit,'String',num2str(delay));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


% --- Executes during object creation, after setting all properties.
function pseudowavelengthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pseudowavelengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pseudofrequencyEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pseudofrequencyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature laser_frequency ...
       wavenumber nyquist spectralResolution ...
       calibration_button fr_real0 f_0 calibration_derivative ...
       apodization_type
% Hints: get(hObject,'String') returns contents of pseudofrequencyEdit as text
%        str2double(get(hObject,'String')) returns contents of pseudofrequencyEdit as a double
   
%update the Edit window
pseudofrequency=str2num(get(handles.pseudofrequencyEdit,'String'));
set(handles.pseudofrequencyEdit,'String',pseudofrequency);

pseudowavelength=pseudofrequency^(-1); %[um]

if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
    %WAVELENGTH AS FUNCTION OF THE PSEUDOWAVELENGTH (OPTICAL CYCLE IN um)
    %wavelength=Delta_n*sin(apexAngle)*pseudowavelength
    %as Delta_n depends on wavelength a numerical method to solve this equation
    %in wavelength is necessary
    %We solve f(wavelength)=0 with the Newton method
    f=@(variable) variable-abs(diff(nCryst(crystal,variable*1e-3,temperature))).*sin(apexAngle/180*pi).*pseudowavelength.*1e3; %[nm]
    initialGuess=Delta_n.*sin(apexAngle/180*pi).*pseudowavelength.*1e3; %initial guess [nm] exploiting the actual Delta_n
    iterations=100; %maximum number of iterations in the method
    tollerance=1e-5; %error tollerance [nm]
    wavelength=newton_method(f,initialGuess,iterations,tollerance); %[nm]

    frequency=c/wavelength*1e-3; %[THz]
    
    Delta_n=abs(diff(nCryst(crystal,wavelength.*1e-3,temperature)));

    delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]
    
else %if calibration is loaded: calculation using experimental TWINS calibration
    frequency=interp1(f_0,fr_real0,pseudofrequency); %[THz]
    wavelength=c/frequency*1e-3; %[nm]
    
    delay=interp1(f_0,calibration_derivative,pseudofrequency).^(-1).*step.*1e3; %[fs] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
end

wavenumber=(laser_frequency-frequency)/c*1e10; %[cm-1]

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.wavelengthEdit,'String',num2str(wavelength));
set(handles.frequencyEdit,'String',num2str(frequency));
set(handles.wavenumberEdit,'String',num2str(wavenumber));
set(handles.pseudowavelengthEdit,'String',num2str(pseudowavelength));
set(handles.delayEdit,'String',num2str(delay));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


% --- Executes during object creation, after setting all properties.
function pseudofrequencyEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pseudofrequencyEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delayEdit_Callback(hObject, eventdata, handles)
% hObject    handle to delayEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature laser_frequency ...
       wavenumber nyquist spectralResolution ...
       calibration_button fr_real0 f_0 calibration_derivative ...
       apodization_type
% Hints: get(hObject,'String') returns contents of delayEdit as text
%        str2double(get(hObject,'String')) returns contents of delayEdit as a double
   
%update the Edit window
delay=str2num(get(handles.delayEdit,'String'));
set(handles.delayEdit,'String',delay);

if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
    
    step=delay*c/(Delta_n*sin(apexAngle/180*pi))*1e-9; %[um]
    
else %if calibration is loaded: calculation using experimental TWINS calibration

    step=interp1(f_0,calibration_derivative,pseudofrequency).*delay.*1e-3; %[um] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX

end

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.stepEdit,'String',num2str(step));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));



% --- Executes during object creation, after setting all properties.
function delayEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delayEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stepEdit_Callback(hObject, eventdata, handles)
% hObject    handle to stepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature laser_frequency ...
       wavenumber nyquist spectralResolution ...
       calibration_button fr_real0 f_0 calibration_derivative ...
       apodization_type
% Hints: get(hObject,'String') returns contents of stepEdit as text
%        str2double(get(hObject,'String')) returns contents of stepEdit as a double
   
%update the Edit window
step=str2num(get(handles.stepEdit,'String'));
set(handles.stepEdit,'String',step);

if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
    
    delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]
    
else %if calibration is loaded: calculation using experimental TWINS calibration
    
    delay=interp1(f_0,calibration_derivative,pseudofrequency).^(-1).*step.*1e3; %[fs] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX

end

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.delayEdit,'String',num2str(delay));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


% --- Executes during object creation, after setting all properties.
function stepEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nyquistEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nyquistEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature laser_frequency ...
       wavenumber nyquist spectralResolution ...
       calibration_button fr_real0 f_0 calibration_derivative ...
       apodization_type
% Hints: get(hObject,'String') returns contents of nyquistEdit as text
%        str2double(get(hObject,'String')) returns contents of nyquistEdit as a double
   
%update the Edit window
nyquist=str2num(get(handles.nyquistEdit,'String'));
set(handles.nyquistEdit,'String',nyquist);

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        step=1/(2*nyquist); %[um]
        if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
            delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]
        else %if calibration is loaded: calculation using experimental TWINS calibration
            delay=interp1(f_0,calibration_derivative,pseudofrequency).^(-1).*step.*1e3; %[fs] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
        end
    case 2 %frequency [THz]
        delay=1/(2*nyquist)*1e3;
        if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
            step=delay/(Delta_n/c*sin(apexAngle/180*pi))*1e-9; %[um]
        else %if calibration is loaded: calculation using experimental TWINS calibration
            step=interp1(f_0,calibration_derivative,pseudofrequency).*delay.*1e-3; %[um] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
        end
    case 3 %wavelength [nm]
        delay=1/(2*(c/nyquist*1e-3))*1e3;
        if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
            step=delay/(Delta_n/c*sin(apexAngle/180*pi))*1e-9; %[um]
        else %if calibration is loaded: calculation using experimental TWINS calibration
            step=interp1(f_0,calibration_derivative,pseudofrequency).*delay.*1e-3; %[um] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
        end
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.delayEdit,'String',num2str(delay));
set(handles.stepEdit,'String',num2str(step));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


% --- Executes during object creation, after setting all properties.
function nyquistEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nyquistEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spResolutionEdit_Callback(hObject, eventdata, handles)
% hObject    handle to spResolutionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature laser_frequency ...
       wavenumber nyquist spectralResolution ...
       calibration_button fr_real0 f_0 calibration_derivative ...
       apodization_type
% Hints: get(hObject,'String') returns contents of spResolutionEdit as text
%        str2double(get(hObject,'String')) returns contents of spResolutionEdit as a double
   
%update the Edit window
spectralResolution=str2num(get(handles.spResolutionEdit,'String'));
set(handles.spResolutionEdit,'String',spectralResolution);

switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        Delta_f=spectralResolution; %[THz]
        
        switch apodization_type
            case 0 %rectangular apodization
                delay=1.2043./Delta_f.*1e3; %[fs] %the conversion factor 1.2043 has been calculated in Matlab
            case 1 %Happ-Genzel apodization
                delay=1.8152./Delta_f.*1e3; %[fs] %the conversion factor 1.8152 has been calculated in Matlab
            case 2 %3-term Blackmann-Harris
                delay=2.2694./Delta_f.*1e3; %[fs] %the conversion factor 2.2694 has been calculated in Matlab
            case 3 %4-term Blackmann-Harris
                delay=2.6611./Delta_f.*1e3; %[fs] %the conversion factor 2.6611 has been calculated in Matlab
            case 4 %Trinagular
                delay=1.7683./Delta_f.*1e3; %[fs] %the conversion factor 1.7683 has been calculated in Matlab
        end

        if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
            step=delay/(Delta_n/c*sin(apexAngle/180*pi))*1e-9; %[um]
        else %if calibration is loaded: calculation using experimental TWINS calibration
            step=interp1(f_0,calibration_derivative,pseudofrequency).*delay.*1e-3; %[um] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
        end
        
    case 2 %wavelength [nm]
        Delta_f=c/wavelength*spectralResolution/(wavelength-spectralResolution)*1e-3; %[THz]
        
        switch apodization_type
            case 0 %rectangular apodization
                delay=1.2043./Delta_f.*1e3; %[fs] %the conversion factor 1.2043 has been calculated in Matlab
            case 1 %Happ-Genzel apodization
                delay=1.8152./Delta_f.*1e3; %[fs] %the conversion factor 1.8152 has been calculated in Matlab
            case 2 %3-term Blackmann-Harris
                delay=2.2694./Delta_f.*1e3; %[fs] %the conversion factor 2.2694 has been calculated in Matlab
            case 3 %4-term Blackmann-Harris
                delay=2.6611./Delta_f.*1e3; %[fs] %the conversion factor 2.6611 has been calculated in Matlab
            case 4 %Trinagular
                delay=1.7683./Delta_f.*1e3; %[fs] %the conversion factor 1.7683 has been calculated in Matlab
        end
        
        if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
            step=delay/(Delta_n/c*sin(apexAngle/180*pi))*1e-9; %[um]
        else %if calibration is loaded: calculation using experimental TWINS calibration
            step=interp1(f_0,calibration_derivative,pseudofrequency).*delay.*1e-3; %[um] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
        end
    case 3 %wavenumber [cm-1]
        Delta_f=spectralResolution*c*1e-10; %[THz]
        
        switch apodization_type
            case 0 %rectangular apodization
                delay=1.2043./Delta_f.*1e3; %[fs] %the conversion factor 1.2043 has been calculated in Matlab
            case 1 %Happ-Genzel apodization
                delay=1.8152./Delta_f.*1e3; %[fs] %the conversion factor 1.8152 has been calculated in Matlab
            case 2 %3-term Blackmann-Harris
                delay=2.2694./Delta_f.*1e3; %[fs] %the conversion factor 2.2694 has been calculated in Matlab
            case 3 %4-term Blackmann-Harris
                delay=2.6611./Delta_f.*1e3; %[fs] %the conversion factor 2.6611 has been calculated in Matlab
            case 4 %Trinagular
                delay=1.7683./Delta_f.*1e3; %[fs] %the conversion factor 1.7683 has been calculated in Matlab
        end
        
        if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
            step=delay/(Delta_n/c*sin(apexAngle/180*pi))*1e-9; %[um]
        else %if calibration is loaded: calculation using experimental TWINS calibration
            step=interp1(f_0,calibration_derivative,pseudofrequency).*delay.*1e-3; %[um] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
        end
    case 4 %pseudofrequency [um^[-1}]
                switch apodization_type
            case 0 %rectangular apodization
                step=1.2043./spectralResolution; %[um] %the conversion factor 1.2043 has been calculated in Matlab
            case 1 %Happ-Genzel apodization
                step=1.8152./spectralResolution; %[um] %the conversion factor 1.8152 has been calculated in Matlab
            case 2 %3-term Blackmann-Harris
                step=2.2694./spectralResolution; %[um] %the conversion factor 2.2694 has been calculated in Matlab
            case 3 %4-term Blackmann-Harris
                step=2.6611./spectralResolution; %[um] %the conversion factor 2.6611 has been calculated in Matlab
            case 4 %Trinagular
                step=1.7683./spectralResolution; %[um] %the conversion factor 1.7683 has been calculated in Matlab
        end

        if calibration_button==0 %if calibration is not loaded: calculation using TWINS theory formula
            delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]
        else %if calibration is loaded: calculation using experimental TWINS calibration
            delay=interp1(f_0,calibration_derivative,pseudofrequency).^(-1).*step.*1e3; %[fs] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
        end
end

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

%update changed parameters
set(handles.delayEdit,'String',num2str(delay));
set(handles.stepEdit,'String',num2str(step));
set(handles.nyquistEdit,'String',num2str(nyquist));



% --- Executes during object creation, after setting all properties.
function spResolutionEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spResolutionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nyquistMenu.
function nyquistMenu_Callback(hObject, eventdata, handles)
% hObject    handle to nyquistMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature nyquist spectralResolution
% Hints: contents = cellstr(get(hObject,'String')) returns nyquistMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nyquistMenu

switch get(hObject,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

%update changed parameters
set(handles.nyquistEdit,'String',num2str(nyquist));


% --- Executes during object creation, after setting all properties.
function nyquistMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nyquistMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in spResolutionMenu.
function spResolutionMenu_Callback(hObject, eventdata, handles)
% hObject    handle to spResolutionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature nyquist spectralResolution ...
       apodization_type
% Hints: contents = cellstr(get(hObject,'String')) returns spResolutionMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spResolutionMenu

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(hObject,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=1/step;
end

%update changed parameters
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


% --- Executes during object creation, after setting all properties.
function spResolutionMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spResolutionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function temperatureEdit_Callback(hObject, eventdata, handles)
% hObject    handle to temperatureEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature nyquist spectralResolution ...
       apodization_type
% Hints: get(hObject,'String') returns contents of temperatureEdit as text
%        str2double(get(hObject,'String')) returns contents of temperatureEdit as a double
   
%update the Edit window
temperature=str2num(get(handles.temperatureEdit,'String'));
set(handles.temperatureEdit,'String',temperature);

Delta_n=abs(diff(nCryst(crystal,wavelength.*1e-3,temperature)));

pseudowavelength=(1/frequency)/(Delta_n/c*sin(apexAngle/180*pi))*1e-6; %optical cycle [um]
pseudofrequency=pseudowavelength^(-1); %[um-1]

delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.pseudowavelengthEdit,'String',num2str(pseudowavelength));
set(handles.pseudofrequencyEdit,'String',num2str(pseudofrequency));
set(handles.delayEdit,'String',num2str(delay));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


% --- Executes during object creation, after setting all properties.
function temperatureEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temperatureEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calibration_enable.
function calibration_enable_Callback(hObject, eventdata, handles)
% hObject    handle to calibration_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global crystal pseudofrequency pseudowavelength wavelength frequency ...
       step delay Delta_n c apexAngle temperature laser_frequency ...
       wavenumber nyquist spectralResolution ...
       calibration_button fr_real0 f_0 calibration_derivative ...
       apodization_type
% Hint: get(hObject,'Value') returns toggle state of calibration_enable

switch get(hObject,'Value')
    case 0 %button false (void)
        calibration_button=0; %no calibration loaded
        set(handles.calibration, 'String', '(none)'); %text in the GUI that indicates loaded calibration
        
        %enable crystal pop-up menu, apex angle box and temperature box
        set(handles.crystalMenu,'Enable','on');
        set(handles.apexAngleEdit,'Enable','on');
        set(handles.temperatureEdit,'Enable','on');
        
        %update parameter according to TWINS theory (no calibration is
        %loaded)
        Delta_n=abs(diff(nCryst(crystal,wavelength.*1e-3,temperature)));
        pseudowavelength=(1/frequency)/(Delta_n/c*sin(apexAngle/180*pi))*1e-6; %optical cycle [um]
        pseudofrequency=pseudowavelength^(-1); %[um-1]
        delay=Delta_n/c*sin(apexAngle/180*pi)*step*1e9; %[fs]
        
    case 1 %button true (filled)
        calibration_button=1; %calibration loaded
        
        [file_calibration,path_calibration] = uigetfile('*.mat', 'Load Calibration'); %load calibration
        filename_calibration=[path_calibration '\' file_calibration];
        
        set(handles.calibration, 'String', file_calibration); %text in the GUI that indicates loaded calibration
        
        input_calibration=load(filename_calibration);
        f_0=input_calibration.f_0; %calibration pseudofrequencies axis [um^-1]
        fr_real0=input_calibration.fr_real0; %calibration frequencies axis [THz]
        
        calibration_derivative=zeros(size(fr_real0));
        calibration_derivative(1:end-1)=diff(fr_real0)./diff(f_0);
        calibration_derivative(end)=calibration_derivative(end-1);
        
        %disable crystal pop-up menu, apex angle box and temperature box
        set(handles.crystalMenu,'Enable','off');
        set(handles.apexAngleEdit,'Enable','off');
        set(handles.temperatureEdit,'Enable','off');
        
        %update parameter according to experimental TWINS calibration (calibration is
        %loaded)
        pseudofrequency=interp1(fr_real0,f_0,frequency); %[um-1]
        pseudowavelength=pseudofrequency^(-1); %[um-1] 
        delay=interp1(f_0,calibration_derivative,pseudofrequency).^(-1).*step.*1e3; %[fs] deltaTau=(d CALIBRATION / d PSEUDOFREQUENCY)^(-1)*deltaX
end

switch get(handles.nyquistMenu,'Value')
    case 1 %pseudofrequency [um^[-1}]
        nyquist=1/(2*step);
    case 2 %frequency [THz]
        nyquist=1/(2*delay)*1e3;
    case 3 %wavelength [nm]
        nyquist=c/(1/(2*delay)*1e3)*1e-3;
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.frequencyEdit,'String',num2str(frequency));
set(handles.wavenumberEdit,'String',num2str(wavenumber));
set(handles.pseudowavelengthEdit,'String',num2str(pseudowavelength));
set(handles.pseudofrequencyEdit,'String',num2str(pseudofrequency));
set(handles.delayEdit,'String',num2str(delay));
set(handles.nyquistEdit,'String',num2str(nyquist));
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


% --- Executes on selection change in apodizationMenu.
function apodizationMenu_Callback(hObject, eventdata, handles)
% hObject    handle to apodizationMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global apodization_type Delta_f c wavelength delay step spectralResolution
% Hints: contents = cellstr(get(hObject,'String')) returns apodizationMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from apodizationMenu

switch get(hObject,'Value')
    case 1
        apodization_type=0; %rectangular
    case 2
        apodization_type=1; %Happ-Genzel
    case 3
        apodization_type=2; %3-term Blackmann-Harris
    case 4
        apodization_type=3; %4-term Blackmann-Harris
    case 5
        apodization_type=4; %Triangular
end

Delta_f=FWHM_apodization(apodization_type,delay).*1e3; %spectral resolution in THz
switch get(handles.spResolutionMenu,'Value')
    case 1 %frequency [THz]
        spectralResolution=Delta_f;
    case 2 %wavelength [nm]
        spectralResolution=(Delta_f/c*1e3)/(1+Delta_f/c*1e3)*wavelength^2;
    case 3 %wavenumber [cm-1]
        spectralResolution=Delta_f/c*1e10;
    case 4 %pseudofrequency [um^[-1}]
        spectralResolution=FWHM_apodization(apodization_type,step);
end

%update changed parameters
set(handles.spResolutionEdit,'String',num2str(spectralResolution));


% --- Executes during object creation, after setting all properties.
function apodizationMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apodizationMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
