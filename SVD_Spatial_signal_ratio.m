function varargout = SVD_Spatial_signal_ratio(varargin)
% SVD_SPATIAL_SIGNAL_RATIO MATLAB code for SVD_Spatial_signal_ratio.fig
%      SVD_SPATIAL_SIGNAL_RATIO, by itself, creates a new SVD_SPATIAL_SIGNAL_RATIO or raises the existing
%      singleton*.
%
%      H = SVD_SPATIAL_SIGNAL_RATIO returns the handle to a new SVD_SPATIAL_SIGNAL_RATIO or the handle to
%      the existing singleton*.
%
%      SVD_SPATIAL_SIGNAL_RATIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SVD_SPATIAL_SIGNAL_RATIO.M with the given input arguments.
%
%      SVD_SPATIAL_SIGNAL_RATIO('Property','Value',...) creates a new SVD_SPATIAL_SIGNAL_RATIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SVD_Spatial_signal_ratio_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SVD_Spatial_signal_ratio_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SVD_Spatial_signal_ratio

% Last Modified by GUIDE v2.5 29-May-2023 15:30:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SVD_Spatial_signal_ratio_OpeningFcn, ...
                   'gui_OutputFcn',  @SVD_Spatial_signal_ratio_OutputFcn, ...
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


% --- Executes just before SVD_Spatial_signal_ratio is made visible.
function SVD_Spatial_signal_ratio_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SVD_Spatial_signal_ratio (see VARARGIN)

% Choose default command line output for SVD_Spatial_signal_ratio
handles.output = hObject;

global U V S a b p_variables spatial_signal_ratio t f threshold x_selected y_selected ...
    apod yLimMin yLimMax xLimMin xLimMax black ...
    saturation gamma max_image immagine spectrum ...
    HyperMatrix_new spectrum_y_fixed spectrum_y_lim interferogram interf_y_lim ...
    interf_y_fixed pixels_mean n_SVs t_min t_max

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SVD_Spatial_signal_ratio wait for user response (see UIRESUME)
% uiwait(handles.figure1);

U=varargin{1}; %U matrix of SVD of HyperMatrix (temporal hypercube)
S=varargin{2}; %S matrix of SVD of HyperMatrix (temporal hypercube)
V=varargin{3}; %V matrix of SVD of HyperMatrix (temporal hypercube)
a=varargin{4}; %this is size(HyperMatrix,1); y dimension
b=varargin{5}; %this is size(HyperMatrix,2); x dimension
p_variables=varargin{6}; %this is size(HyperMatrix,3); delay dimension
spatial_signal_ratio=varargin{7}; %spatial signal ratio vector
threshold=varargin{8}; %first threshold for Spatial Signal Ratio to be considered
n_SVs=varargin{9}; %number of Singular Values (SVs) considered
t=varargin{10}; %delay (or motor position axis)
f=varargin{11}; %frequency (or pseudofrequency axis)

if threshold>max(spatial_signal_ratio)
    threshold=max(spatial_signal_ratio);
end

if threshold<min(spatial_signal_ratio)
    threshold=min(spatial_signal_ratio);
end

%range for y axis in Spatial Signal Ratio Plot
yLimMin=min(spatial_signal_ratio).*0.9;
yLimMax=max(spatial_signal_ratio).*1.1;

%range for x axis in Spatial Signal Ratio Plot
xLimMin=0;
xLimMax=n_SVs;

set(handles.thresholdValue,'String',num2str(threshold));
set(handles.YLimMinValue,'String',num2str(yLimMin));
set(handles.YLimMaxValue,'String',num2str(yLimMax));
set(handles.N_SVs_Edit,'String',num2str(n_SVs));

plot(handles.Plot_Spatial_signal_ratio,spatial_signal_ratio(1:n_SVs),'b','linewidth',2); axis tight;
set(handles.Plot_Spatial_signal_ratio,'YLim',[yLimMin yLimMax]);
set(handles.Plot_Spatial_signal_ratio,'XLim',[xLimMin xLimMax]);

hold on;
yline(handles.Plot_Spatial_signal_ratio,threshold,'r','linewidth',1);
hold off; %this is important to manage two different graphs in the GUI

SVs=diag(S); %singular values (SVs)

SVs_swap=zeros(size(SVs));
SVs_swap(1:end-(length(SVs)-n_SVs))=SVs(1:n_SVs); %keep only the first n_SVs SVs

SVs=SVs_swap;
clear SVs_swap

SVs(spatial_signal_ratio<threshold)=0; %put to 0 all the SVs associated to low spatial_signal_ratio

S_new=diag(SVs); %recalculate the S matrix

X_new=(U*S_new*V')'; %recalculate the denoised dataset

HyperMatrix_new=reshape(X_new,a,b,p_variables); %reshape properly the dataset

Sample=sum(squeeze(abs(HyperMatrix_new)),3);
immagine=zeros(size(Sample,1),size(Sample,2),3);
immagine(:,:,1)=Sample;
immagine(:,:,2)=Sample;
immagine(:,:,3)=Sample;
black=0;
saturation=1;
gamma=1;
max_image=max(max(max(immagine)));

set(handles.SaturationValue,'String',num2str(saturation));
set(handles.BlackValue,'String',num2str(black));
set(handles.GammaValue,'String',num2str(gamma));

immagine_plot=(abs(immagine./max_image-black)/saturation).^gamma;

x_selected=round(size(immagine_plot,2)./2);
y_selected=x_selected;

immagine_plot(y_selected,1:end,1)=0;
immagine_plot(y_selected,1:end,2)=0.7;
immagine_plot(y_selected,1:end,3)=0;
immagine_plot(:,x_selected,1)=0;
immagine_plot(:,x_selected,2)=0.7;
immagine_plot(:,x_selected,3)=0;

image(immagine_plot,'Parent',handles.Immagine); axis equal; axis off;

set(handles.x_selectedValue,'String',num2str(x_selected));
set(handles.y_selectedValue,'String',num2str(y_selected));

pixels_mean=1;
set(handles.NpixelsMeanValue,'String',num2str(pixels_mean));

%interferogram=squeeze(HyperMatrix_new(y_selected,x_selected,:));

if mod(pixels_mean,2)==0 %pixels_mean even
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-(pixels_mean/2-1):y_selected+pixels_mean/2,x_selected-(pixels_mean/2-1):x_selected+pixels_mean/2,:),[1,2]))';
else %pixels_mean odd
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-floor(pixels_mean/2):y_selected+floor(pixels_mean/2),x_selected-floor(pixels_mean/2):x_selected+floor(pixels_mean/2),:),[1,2]))';
end

interferogram=interferogram-mean(interferogram);

t_min=min(t); %minimum of delay t axis considered for FT
t_max=max(t); %maximum of delay t axis considered for FT

set(handles.t_minValue,'String',t_min);
set(handles.t_maxValue,'String',t_max);
set(handles.min_t_value, 'String', num2str(round(t_min,1))); %min t value indicated in the GUI
set(handles.max_t_value, 'String', num2str(round(t_max,1))); %max t value indicated in the GUI

[~,pos_min]=min(abs(t-t_min));
[~,pos_max]=min(abs(t-t_max));

apod=Apodization(1,length(t(pos_min:pos_max)),length(t(pos_min:pos_max))./2,7); %interferogram apodization

%spectrum=abs(FourierDir(t,interferogram'.*apod,f));
spectrum=abs(FourierDir(t(pos_min:pos_max),interferogram(pos_min:pos_max).*apod,f));

spectrum_y_fixed=0; %condition 'do not fix spectrum y-axis'
set(handles.FixYSpectrum,'Value',spectrum_y_fixed);

interf_y_fixed=0; %condition 'do not fix interferogram y-axis'
set(handles.FixYInterferogram,'Value',interf_y_fixed);

plot(handles.Plot_spectr,f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_spectr,'XLim',[min(f) max(f)]);
spectrum_y_lim=get(handles.Plot_spectr,'YLim');

plot(handles.Plot_Interf,t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_Interf,'XLim',[t_min t_max]);
interf_y_lim=get(handles.Plot_Interf,'YLim');





% --- Outputs from this function are returned to the command line.
function varargout = SVD_Spatial_signal_ratio_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function thresholdValue_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdValue as text
%        str2double(get(hObject,'String')) returns contents of thresholdValue as a double

global threshold yLimMin yLimMax xLimMin xLimMax ...
    spatial_signal_ratio U S V a b p_variables ...
    black saturation gamma max_image x_selected y_selected t f apod immagine ...
    spectrum HyperMatrix_new spectrum_y_lim spectrum_y_fixed interferogram ...
    interf_y_fixed interf_y_lim pixels_mean n_SVs t_min t_max

threshold=str2num(get(handles.thresholdValue,'String'));

if threshold>max(spatial_signal_ratio)
    threshold=max(spatial_signal_ratio);
end

if threshold<min(spatial_signal_ratio)
    threshold=min(spatial_signal_ratio);
end

set(handles.thresholdValue,'String',threshold);

hhh=waitbar(1/2,'Please, wait...');

plot(handles.Plot_Spatial_signal_ratio,spatial_signal_ratio(1:n_SVs),'b','linewidth',2); axis tight;
set(handles.Plot_Spatial_signal_ratio,'YLim',[yLimMin yLimMax]);
set(handles.Plot_Spatial_signal_ratio,'XLim',[xLimMin xLimMax]);

hold on;
yline(handles.Plot_Spatial_signal_ratio,threshold,'r','linewidth',1);
hold off; %this is important to manage two different graphs in the GUI

SVs=diag(S); %singular values (SVs)

SVs_swap=zeros(size(SVs));
SVs_swap(1:end-(length(SVs)-n_SVs))=SVs(1:n_SVs); %keep only the first n_SVs SVs

SVs=SVs_swap;
clear SVs_swap

SVs(spatial_signal_ratio<threshold)=0; %put to 0 all the SVs associated to low spatial_signal_ratio

S_new=diag(SVs); %recalculate the S matrix

X_new=(U*S_new*V')'; %recalculate the denoised dataset

HyperMatrix_new=reshape(X_new,a,b,p_variables); %reshape properly the dataset

waitbar(1,hhh,'Please, wait...');
close(hhh);

Sample=sum(squeeze(abs(HyperMatrix_new)),3);
immagine=zeros(size(Sample,1),size(Sample,2),3);
immagine(:,:,1)=Sample;
immagine(:,:,2)=Sample;
immagine(:,:,3)=Sample;

max_image=max(max(max(immagine)));

immagine_plot=(abs(immagine./max_image-black)/saturation).^gamma;

immagine_plot(y_selected,1:end,1)=0;
immagine_plot(y_selected,1:end,2)=0.7;
immagine_plot(y_selected,1:end,3)=0;
immagine_plot(:,x_selected,1)=0;
immagine_plot(:,x_selected,2)=0.7;
immagine_plot(:,x_selected,3)=0;

image(immagine_plot,'Parent',handles.Immagine); axis equal; axis off;

%interferogram=squeeze(HyperMatrix_new(y_selected,x_selected,:));
if mod(pixels_mean,2)==0 %pixels_mean even
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-(pixels_mean/2-1):y_selected+pixels_mean/2,x_selected-(pixels_mean/2-1):x_selected+pixels_mean/2,:),[1,2]))';
else %pixels_mean odd
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-floor(pixels_mean/2):y_selected+floor(pixels_mean/2),x_selected-floor(pixels_mean/2):x_selected+floor(pixels_mean/2),:),[1,2]))';
end
interferogram=interferogram-mean(interferogram);

[~,pos_min]=min(abs(t-t_min));
[~,pos_max]=min(abs(t-t_max));

apod=Apodization(1,length(t(pos_min:pos_max)),length(t(pos_min:pos_max))./2,7); %interferogram apodization
%spectrum=abs(FourierDir(t,interferogram'.*apod,f));
spectrum=abs(FourierDir(t(pos_min:pos_max),interferogram(pos_min:pos_max).*apod,f));


plot(handles.Plot_spectr,f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_spectr,'XLim',[min(f) max(f)]);

if spectrum_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_spectr,'YLim',spectrum_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    spectrum_y_lim=get(handles.Plot_spectr,'YLim'); %take the new limits
end

plot(handles.Plot_Interf,t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_Interf,'XLim',[t_min t_max]);

if interf_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_Interf,'YLim',interf_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    interf_y_lim=get(handles.Plot_Interf,'YLim'); %take the new limits
end






% --- Executes during object creation, after setting all properties.
function thresholdValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_selectedValue_Callback(hObject, eventdata, handles)
% hObject    handle to x_selectedValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_selectedValue as text
%        str2double(get(hObject,'String')) returns contents of x_selectedValue as a double

global saturation black gamma immagine max_image x_selected y_selected ...
    spectrum f t apod HyperMatrix_new spectrum_y_fixed spectrum_y_lim ...
    interferogram interf_y_fixed interf_y_lim pixels_mean t_min t_max

x_selected=str2num(get(handles.x_selectedValue,'String'));

if x_selected<1
    x_selected=1;
end

if x_selected>size(immagine,2)
    x_selected=size(immagine,2);
end

set(handles.x_selectedValue,'String',x_selected);

if y_selected+floor(pixels_mean./2)>size(immagine,1)
    pixels_mean=2.*(size(immagine,1)-y_selected)+1; %the pixel mean square cannot overcome the lower image limit
end

if y_selected-floor(pixels_mean./2)<1
    pixels_mean=2.*y_selected-1; %the pixel mean square cannot overcome the upper image limit
end

if x_selected+floor(pixels_mean./2)>size(immagine,2)
    pixels_mean=2.*(size(immagine,2)-x_selected)+1; %the pixel mean square cannot overcome the right image limit
end

if x_selected-floor(pixels_mean./2)<1
    pixels_mean=2.*x_selected-1; %the pixel mean square cannot overcome the left image limit
end

set(handles.NpixelsMeanValue,'String',pixels_mean);

immagine_plot=(abs(immagine./max_image-black)/saturation).^gamma;

immagine_plot(y_selected,1:end,1)=0;
immagine_plot(y_selected,1:end,2)=0.7;
immagine_plot(y_selected,1:end,3)=0;
immagine_plot(:,x_selected,1)=0;
immagine_plot(:,x_selected,2)=0.7;
immagine_plot(:,x_selected,3)=0;

image(immagine_plot,'Parent',handles.Immagine); axis equal; axis off;

%interferogram=squeeze(HyperMatrix_new(y_selected,x_selected,:));
if mod(pixels_mean,2)==0 %pixels_mean even
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-(pixels_mean/2-1):y_selected+pixels_mean/2,x_selected-(pixels_mean/2-1):x_selected+pixels_mean/2,:),[1,2]))';
else %pixels_mean odd
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-floor(pixels_mean/2):y_selected+floor(pixels_mean/2),x_selected-floor(pixels_mean/2):x_selected+floor(pixels_mean/2),:),[1,2]))';
end
interferogram=interferogram-mean(interferogram);

[~,pos_min]=min(abs(t-t_min));
[~,pos_max]=min(abs(t-t_max));

apod=Apodization(1,length(t(pos_min:pos_max)),length(t(pos_min:pos_max))./2,7); %interferogram apodization

%spectrum=abs(FourierDir(t,interferogram'.*apod,f));
spectrum=abs(FourierDir(t(pos_min:pos_max),interferogram(pos_min:pos_max).*apod,f));

plot(handles.Plot_spectr,f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_spectr,'XLim',[min(f) max(f)]);

if spectrum_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_spectr,'YLim',spectrum_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    spectrum_y_lim=get(handles.Plot_spectr,'YLim'); %take the new limits
end

plot(handles.Plot_Interf,t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_Interf,'XLim',[t_min t_max]);

if interf_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_Interf,'YLim',interf_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    interf_y_lim=get(handles.Plot_Interf,'YLim'); %take the new limits
end

% --- Executes during object creation, after setting all properties.
function x_selectedValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_selectedValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_selectedValue_Callback(hObject, eventdata, handles)
% hObject    handle to y_selectedValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_selectedValue as text
%        str2double(get(hObject,'String')) returns contents of y_selectedValue as a double

global saturation black gamma immagine max_image x_selected y_selected ...
    spectrum f t apod HyperMatrix_new spectrum_y_fixed spectrum_y_lim ...
    interferogram interf_y_fixed interf_y_lim pixels_mean t_min t_max

y_selected=str2num(get(handles.y_selectedValue,'String'));

if y_selected<1
    y_selected=1;
end

if y_selected>size(immagine,2)
    y_selected=size(immagine,2);
end

set(handles.y_selectedValue,'String',y_selected);

if y_selected+floor(pixels_mean./2)>size(immagine,1)
    pixels_mean=2.*(size(immagine,1)-y_selected)+1; %the pixel mean square cannot overcome the lower image limit
end

if y_selected-floor(pixels_mean./2)<1
    pixels_mean=2.*y_selected-1; %the pixel mean square cannot overcome the upper image limit
end

if x_selected+floor(pixels_mean./2)>size(immagine,2)
    pixels_mean=2.*(size(immagine,2)-x_selected)+1; %the pixel mean square cannot overcome the right image limit
end

if x_selected-floor(pixels_mean./2)<1
    pixels_mean=2.*x_selected-1; %the pixel mean square cannot overcome the left image limit
end

set(handles.NpixelsMeanValue,'String',pixels_mean);

immagine_plot=(abs(immagine./max_image-black)/saturation).^gamma;

immagine_plot(y_selected,1:end,1)=0;
immagine_plot(y_selected,1:end,2)=0.7;
immagine_plot(y_selected,1:end,3)=0;
immagine_plot(:,x_selected,1)=0;
immagine_plot(:,x_selected,2)=0.7;
immagine_plot(:,x_selected,3)=0;

image(immagine_plot,'Parent',handles.Immagine); axis equal; axis off;

%interferogram=squeeze(HyperMatrix_new(y_selected,x_selected,:));
if mod(pixels_mean,2)==0 %pixels_mean even
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-(pixels_mean/2-1):y_selected+pixels_mean/2,x_selected-(pixels_mean/2-1):x_selected+pixels_mean/2,:),[1,2]))';
else %pixels_mean odd
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-floor(pixels_mean/2):y_selected+floor(pixels_mean/2),x_selected-floor(pixels_mean/2):x_selected+floor(pixels_mean/2),:),[1,2]))';
end
interferogram=interferogram-mean(interferogram);

[~,pos_min]=min(abs(t-t_min));
[~,pos_max]=min(abs(t-t_max));

apod=Apodization(1,length(t(pos_min:pos_max)),length(t(pos_min:pos_max))./2,7); %interferogram apodization

%spectrum=abs(FourierDir(t,interferogram'.*apod,f));
spectrum=abs(FourierDir(t(pos_min:pos_max),interferogram(pos_min:pos_max).*apod,f));

plot(handles.Plot_spectr,f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_spectr,'XLim',[min(f) max(f)]);

if spectrum_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_spectr,'YLim',spectrum_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    spectrum_y_lim=get(handles.Plot_spectr,'YLim'); %take the new limits
end

plot(handles.Plot_Interf,t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_Interf,'XLim',[t_min t_max]);

if interf_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_Interf,'YLim',interf_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    interf_y_lim=get(handles.Plot_Interf,'YLim'); %take the new limits
end

% --- Executes during object creation, after setting all properties.
function y_selectedValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_selectedValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YLimMinValue_Callback(hObject, eventdata, handles)
% hObject    handle to YLimMinValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YLimMinValue as text
%        str2double(get(hObject,'String')) returns contents of YLimMinValue as a double

global yLimMin yLimMax

yLimMin=str2num(get(handles.YLimMinValue,'String'));

if yLimMin>=yLimMax
    yLimMin=yLimMax.*0.9;
end

set(handles.YLimMinValue,'String',yLimMin);

set(handles.Plot_Spatial_signal_ratio,'YLim',[yLimMin yLimMax]);




% --- Executes during object creation, after setting all properties.
function YLimMinValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YLimMinValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YLimMaxValue_Callback(hObject, eventdata, handles)
% hObject    handle to YLimMaxValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YLimMaxValue as text
%        str2double(get(hObject,'String')) returns contents of YLimMaxValue as a double

global yLimMin yLimMax

yLimMax=str2num(get(handles.YLimMaxValue,'String'));

if yLimMax<=yLimMin
    yLimMax=yLimMin.*1.1;
end

set(handles.YLimMaxValue,'String',yLimMax);

set(handles.Plot_Spatial_signal_ratio,'YLim',[yLimMin yLimMax]);


% --- Executes during object creation, after setting all properties.
function YLimMaxValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YLimMaxValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SaturationValue_Callback(hObject, eventdata, handles)
% hObject    handle to SaturationValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SaturationValue as text
%        str2double(get(hObject,'String')) returns contents of SaturationValue as a double

global saturation black gamma immagine max_image x_selected y_selected ...
    spectrum f spectrum_y_fixed spectrum_y_lim t interferogram interf_y_fixed ...
    interf_y_lim t_min t_max

saturation=str2num(get(handles.SaturationValue,'String'));

if saturation<=0
    saturation=0.1;
end

if saturation>1
    saturation=1;
end

set(handles.SaturationValue,'String',saturation);

immagine_plot=(abs(immagine./max_image-black)/saturation).^gamma;

immagine_plot(y_selected,1:end,1)=0;
immagine_plot(y_selected,1:end,2)=0.7;
immagine_plot(y_selected,1:end,3)=0;
immagine_plot(:,x_selected,1)=0;
immagine_plot(:,x_selected,2)=0.7;
immagine_plot(:,x_selected,3)=0;

image(immagine_plot,'Parent',handles.Immagine); axis equal; axis off;

%plot again spectrum otherwise Matlab shrink it
plot(handles.Plot_spectr,f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_spectr,'XLim',[min(f) max(f)]);

if spectrum_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_spectr,'YLim',spectrum_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    spectrum_y_lim=get(handles.Plot_spectr,'YLim'); %take the new limits
end

plot(handles.Plot_Interf,t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_Interf,'XLim',[t_min t_max]);

if interf_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_Interf,'YLim',interf_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    interf_y_lim=get(handles.Plot_Interf,'YLim'); %take the new limits
end



% --- Executes during object creation, after setting all properties.
function SaturationValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SaturationValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BlackValue_Callback(hObject, eventdata, handles)
% hObject    handle to BlackValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BlackValue as text
%        str2double(get(hObject,'String')) returns contents of BlackValue as a double

global saturation black gamma immagine max_image x_selected y_selected ...
    spectrum f spectrum_y_fixed spectrum_y_lim t interferogram interf_y_fixed ...
    interf_y_lim t_min t_max

black=str2num(get(handles.BlackValue,'String'));

if black<0
    black=0;
end

if black>1
    black=1;
end

set(handles.BlackValue,'String',black);

immagine_plot=(abs(immagine./max_image-black)/saturation).^gamma;

immagine_plot(y_selected,1:end,1)=0;
immagine_plot(y_selected,1:end,2)=0.7;
immagine_plot(y_selected,1:end,3)=0;
immagine_plot(:,x_selected,1)=0;
immagine_plot(:,x_selected,2)=0.7;
immagine_plot(:,x_selected,3)=0;

image(immagine_plot,'Parent',handles.Immagine); axis equal; axis off;

%plot again spectrum otherwise Matlab shrink it
plot(handles.Plot_spectr,f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_spectr,'XLim',[min(f) max(f)]);

if spectrum_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_spectr,'YLim',spectrum_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    spectrum_y_lim=get(handles.Plot_spectr,'YLim'); %take the new limits
end

plot(handles.Plot_Interf,t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_Interf,'XLim',[t_min t_max]);

if interf_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_Interf,'YLim',interf_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    interf_y_lim=get(handles.Plot_Interf,'YLim'); %take the new limits
end

% --- Executes during object creation, after setting all properties.
function BlackValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BlackValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function GammaValue_Callback(hObject, eventdata, handles)
% hObject    handle to GammaValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GammaValue as text
%        str2double(get(hObject,'String')) returns contents of GammaValue as a double

global saturation black gamma immagine max_image x_selected y_selected ...
    spectrum f spectrum_y_fixed spectrum_y_lim t interferogram interf_y_fixed ...
    interf_y_lim t_min t_max

gamma=str2num(get(handles.GammaValue,'String'));

if gamma<0
    gamma=0;
end

if gamma>1
    gamma=1;
end

set(handles.GammaValue,'String',gamma);

immagine_plot=(abs(immagine./max_image-black)/saturation).^gamma;

immagine_plot(y_selected,1:end,1)=0;
immagine_plot(y_selected,1:end,2)=0.7;
immagine_plot(y_selected,1:end,3)=0;
immagine_plot(:,x_selected,1)=0;
immagine_plot(:,x_selected,2)=0.7;
immagine_plot(:,x_selected,3)=0;

image(immagine_plot,'Parent',handles.Immagine); axis equal; axis off;

%plot again spectrum otherwise Matlab shrink it
plot(handles.Plot_spectr,f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_spectr,'XLim',[min(f) max(f)]);

if spectrum_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_spectr,'YLim',spectrum_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    spectrum_y_lim=get(handles.Plot_spectr,'YLim'); %take the new limits
end

plot(handles.Plot_Interf,t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_Interf,'XLim',[t_min t_max]);

if interf_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_Interf,'YLim',interf_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    interf_y_lim=get(handles.Plot_Interf,'YLim'); %take the new limits
end

% --- Executes during object creation, after setting all properties.
function GammaValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GammaValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SpectrumInspection.
function SpectrumInspection_Callback(hObject, eventdata, handles)
% hObject    handle to SpectrumInspection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global spectrum f spectrum_y_lim threshold x_selected y_selected pixels_mean ...
    n_SVs

figure;
plot(f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
xlim([min(f) max(f)]);
ylim(spectrum_y_lim);
xlabel('Pseudofrequency [\mum^{-1}]');
title(['Pixel (x=' num2str(x_selected-floor(pixels_mean/2)) ':' num2str(x_selected+floor(pixels_mean/2)) '; y=' num2str(y_selected-floor(pixels_mean/2)) ':' num2str(y_selected+floor(pixels_mean/2)) '); Threshold=' num2str(threshold) '; N SVs=' num2str(n_SVs)]);

WConversion_cryst;


% --- Executes on button press in FixYSpectrum.
function FixYSpectrum_Callback(hObject, eventdata, handles)
% hObject    handle to FixYSpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FixYSpectrum

global spectrum_y_fixed

switch get(hObject,'Value')
    case 0 %button false (void)
        spectrum_y_fixed=0;
    case 1 %button true (filled)
        spectrum_y_fixed=1;
end


% --- Executes on button press in ExitButton.
function ExitButton_Callback(hObject, eventdata, handles)
% hObject    handle to ExitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global threshold
%assignin('caller','threshold',threshold);
%assignin('base','threshold',threshold); %assign to Matlab workspace
close('SVD_Spatial_signal_ratio') %close the GUI



function NpixelsMeanValue_Callback(hObject, eventdata, handles)
% hObject    handle to NpixelsMeanValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NpixelsMeanValue as text
%        str2double(get(hObject,'String')) returns contents of NpixelsMeanValue as a double
global immagine x_selected y_selected ...
    spectrum f t apod HyperMatrix_new spectrum_y_fixed spectrum_y_lim ...
    interferogram interf_y_fixed interf_y_lim pixels_mean t_min t_max

pixels_mean=str2num(get(handles.NpixelsMeanValue,'String'));

if pixels_mean<1
    pixels_mean=1;
end

if y_selected+floor(pixels_mean./2)>size(immagine,1)
    pixels_mean=2.*(size(immagine,1)-y_selected)+1; %the pixel mean square cannot overcome the lower image limit
end

if y_selected-floor(pixels_mean./2)<1
    pixels_mean=2.*y_selected-1; %the pixel mean square cannot overcome the upper image limit
end

if x_selected+floor(pixels_mean./2)>size(immagine,2)
    pixels_mean=2.*(size(immagine,2)-x_selected)+1; %the pixel mean square cannot overcome the right image limit
end

if x_selected-floor(pixels_mean./2)<1
    pixels_mean=2.*x_selected-1; %the pixel mean square cannot overcome the left image limit
end

set(handles.NpixelsMeanValue,'String',pixels_mean);

%interferogram=squeeze(HyperMatrix_new(y_selected,x_selected,:));
if mod(pixels_mean,2)==0 %pixels_mean even
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-(pixels_mean/2-1):y_selected+pixels_mean/2,x_selected-(pixels_mean/2-1):x_selected+pixels_mean/2,:),[1,2]))';
else %pixels_mean odd
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-floor(pixels_mean/2):y_selected+floor(pixels_mean/2),x_selected-floor(pixels_mean/2):x_selected+floor(pixels_mean/2),:),[1,2]))';
end

interferogram=interferogram-mean(interferogram);

[~,pos_min]=min(abs(t-t_min));
[~,pos_max]=min(abs(t-t_max));

apod=Apodization(1,length(t(pos_min:pos_max)),length(t(pos_min:pos_max))./2,7); %interferogram apodization

%spectrum=abs(FourierDir(t,interferogram'.*apod,f));
spectrum=abs(FourierDir(t(pos_min:pos_max),interferogram(pos_min:pos_max).*apod,f));

plot(handles.Plot_spectr,f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_spectr,'XLim',[min(f) max(f)]);

if spectrum_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_spectr,'YLim',spectrum_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    spectrum_y_lim=get(handles.Plot_spectr,'YLim'); %take the new limits
end

plot(handles.Plot_Interf,t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_Interf,'XLim',[t_min t_max]);

if interf_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_Interf,'YLim',interf_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    interf_y_lim=get(handles.Plot_Interf,'YLim'); %take the new limits
end

% --- Executes during object creation, after setting all properties.
function NpixelsMeanValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NpixelsMeanValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in InterferogramInspection.
function InterferogramInspection_Callback(hObject, eventdata, handles)
% hObject    handle to InterferogramInspection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global interferogram t interf_y_lim threshold x_selected y_selected pixels_mean ...
    n_SVs t_min t_max

figure;
plot(t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
xlim([t_min t_max]);
ylim(interf_y_lim);
xlabel('Pseudodelay [\mum]');
title(['Pixel (x=' num2str(x_selected-floor(pixels_mean/2)) ':' num2str(x_selected+floor(pixels_mean/2)) '; y=' num2str(y_selected-floor(pixels_mean/2)) ':' num2str(y_selected+floor(pixels_mean/2)) '); Threshold=' num2str(threshold) '; N SVs=' num2str(n_SVs)]);


% --- Executes on button press in FixYInterferogram.
function FixYInterferogram_Callback(hObject, eventdata, handles)
% hObject    handle to FixYInterferogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FixYInterferogram
global interf_y_fixed

switch get(hObject,'Value')
    case 0 %button false (void)
        interf_y_fixed=0;
    case 1 %button true (filled)
        interf_y_fixed=1;
end



function N_SVs_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to N_SVs_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_SVs_Edit as text
%        str2double(get(hObject,'String')) returns contents of N_SVs_Edit as a double
global threshold yLimMin yLimMax xLimMin xLimMax ...
    spatial_signal_ratio U S V a b p_variables ...
    black saturation gamma max_image x_selected y_selected t f apod immagine ...
    spectrum HyperMatrix_new spectrum_y_lim spectrum_y_fixed interferogram ...
    interf_y_fixed interf_y_lim pixels_mean n_SVs t_min t_max

n_SVs=str2num(get(handles.N_SVs_Edit,'String'));

if n_SVs>length(diag(S)) %number of SVs cannot be larger than the size of S matrix
    n_SVs=length(diag(S));
end

if n_SVs<1
    n_SVs=1;
end

set(handles.N_SVs_Edit,'String',n_SVs);

hhh=waitbar(1/2,'Please, wait...');

plot(handles.Plot_Spatial_signal_ratio,spatial_signal_ratio(1:n_SVs),'b','linewidth',2); axis tight;
set(handles.Plot_Spatial_signal_ratio,'YLim',[yLimMin yLimMax]);
set(handles.Plot_Spatial_signal_ratio,'XLim',[xLimMin xLimMax]);

hold on;
yline(handles.Plot_Spatial_signal_ratio,threshold,'r','linewidth',1);
hold off; %this is important to manage two different graphs in the GUI

SVs=diag(S); %singular values (SVs)

SVs_swap=zeros(size(SVs));
SVs_swap(1:end-(length(SVs)-n_SVs))=SVs(1:n_SVs); %keep only the first n_SVs SVs

SVs=SVs_swap;
clear SVs_swap

SVs(spatial_signal_ratio<threshold)=0; %put to 0 all the SVs associated to low spatial_signal_ratio

S_new=diag(SVs); %recalculate the S matrix

X_new=(U*S_new*V')'; %recalculate the denoised dataset

HyperMatrix_new=reshape(X_new,a,b,p_variables); %reshape properly the dataset

waitbar(1,hhh,'Please, wait...');
close(hhh);

Sample=sum(squeeze(abs(HyperMatrix_new)),3);
immagine=zeros(size(Sample,1),size(Sample,2),3);
immagine(:,:,1)=Sample;
immagine(:,:,2)=Sample;
immagine(:,:,3)=Sample;

max_image=max(max(max(immagine)));

immagine_plot=(abs(immagine./max_image-black)/saturation).^gamma;

immagine_plot(y_selected,1:end,1)=0;
immagine_plot(y_selected,1:end,2)=0.7;
immagine_plot(y_selected,1:end,3)=0;
immagine_plot(:,x_selected,1)=0;
immagine_plot(:,x_selected,2)=0.7;
immagine_plot(:,x_selected,3)=0;

image(immagine_plot,'Parent',handles.Immagine); axis equal; axis off;

%interferogram=squeeze(HyperMatrix_new(y_selected,x_selected,:));
if mod(pixels_mean,2)==0 %pixels_mean even
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-(pixels_mean/2-1):y_selected+pixels_mean/2,x_selected-(pixels_mean/2-1):x_selected+pixels_mean/2,:),[1,2]))';
else %pixels_mean odd
    interferogram=squeeze(mean(HyperMatrix_new(y_selected-floor(pixels_mean/2):y_selected+floor(pixels_mean/2),x_selected-floor(pixels_mean/2):x_selected+floor(pixels_mean/2),:),[1,2]))';
end
interferogram=interferogram-mean(interferogram);

[~,pos_min]=min(abs(t-t_min));
[~,pos_max]=min(abs(t-t_max));

apod=Apodization(1,length(t(pos_min:pos_max)),length(t(pos_min:pos_max))./2,7); %interferogram apodization

%spectrum=abs(FourierDir(t,interferogram'.*apod,f));
spectrum=abs(FourierDir(t(pos_min:pos_max),interferogram(pos_min:pos_max).*apod,f));


plot(handles.Plot_spectr,f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_spectr,'XLim',[min(f) max(f)]);

if spectrum_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_spectr,'YLim',spectrum_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    spectrum_y_lim=get(handles.Plot_spectr,'YLim'); %take the new limits
end

plot(handles.Plot_Interf,t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_Interf,'XLim',[t_min t_max]);

if interf_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_Interf,'YLim',interf_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    interf_y_lim=get(handles.Plot_Interf,'YLim'); %take the new limits
end

% --- Executes during object creation, after setting all properties.
function N_SVs_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_SVs_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function XLimMinValue_Callback(hObject, eventdata, handles)
% hObject    handle to XLimMinValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XLimMinValue as text
%        str2double(get(hObject,'String')) returns contents of XLimMinValue as a double
global xLimMin xLimMax

xLimMin=str2num(get(handles.XLimMinValue,'String'));

if xLimMin>=xLimMax
    xLimMin=xLimMax.*0.9;
end

set(handles.XLimMinValue,'String',xLimMin);

set(handles.Plot_Spatial_signal_ratio,'XLim',[xLimMin xLimMax]);

% --- Executes during object creation, after setting all properties.
function XLimMinValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XLimMinValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function XLimMaxValue_Callback(hObject, eventdata, handles)
% hObject    handle to XLimMaxValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XLimMaxValue as text
%        str2double(get(hObject,'String')) returns contents of XLimMaxValue as a double
global xLimMin xLimMax

xLimMax=str2num(get(handles.XLimMaxValue,'String'));

if xLimMax<=xLimMin
    xLimMax=xLimMin.*1.1;
end

set(handles.XLimMinValue,'String',xLimMin);

set(handles.Plot_Spatial_signal_ratio,'XLim',[xLimMin xLimMax]);

% --- Executes during object creation, after setting all properties.
function XLimMaxValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XLimMaxValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_minValue_Callback(hObject, eventdata, handles)
% hObject    handle to t_minValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_minValue as text
%        str2double(get(hObject,'String')) returns contents of t_minValue as a double
global t_min t_max interferogram t apod spectrum_y_lim spectrum_y_fixed ...
    f spectrum interf_y_fixed interf_y_lim

t_min=str2num(get(handles.t_minValue,'String'));

[~,pos_min]=min(abs(t-t_min));
[~,pos_max]=min(abs(t-t_max));

t_min=t(pos_min);

set(handles.t_minValue,'String',t_min);

apod=Apodization(1,length(t(pos_min:pos_max)),length(t(pos_min:pos_max))./2,7); %interferogram apodization

%spectrum=abs(FourierDir(t,interferogram'.*apod,f));
spectrum=abs(FourierDir(t(pos_min:pos_max),interferogram(pos_min:pos_max).*apod,f));

spectrum_y_fixed=0; %condition 'do not fix spectrum y-axis'
set(handles.FixYSpectrum,'Value',spectrum_y_fixed);

interf_y_fixed=0; %condition 'do not fix interferogram y-axis'
set(handles.FixYInterferogram,'Value',interf_y_fixed);

plot(handles.Plot_spectr,f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_spectr,'XLim',[min(f) max(f)]);

if spectrum_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_spectr,'YLim',spectrum_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    spectrum_y_lim=get(handles.Plot_spectr,'YLim'); %take the new limits
end

plot(handles.Plot_Interf,t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_Interf,'XLim',[t_min t_max]);

if interf_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_Interf,'YLim',interf_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    interf_y_lim=get(handles.Plot_Interf,'YLim'); %take the new limits
end

% --- Executes during object creation, after setting all properties.
function t_minValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_minValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_maxValue_Callback(hObject, eventdata, handles)
% hObject    handle to t_maxValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_maxValue as text
%        str2double(get(hObject,'String')) returns contents of t_maxValue as a double
global t_min t_max interferogram t apod spectrum_y_lim spectrum_y_fixed ...
    f spectrum interf_y_fixed interf_y_lim

t_max=str2num(get(handles.t_maxValue,'String'));

[~,pos_min]=min(abs(t-t_min));
[~,pos_max]=min(abs(t-t_max));

t_max=t(pos_max);

set(handles.t_maxValue,'String',t_max);

apod=Apodization(1,length(t(pos_min:pos_max)),length(t(pos_min:pos_max))./2,7); %interferogram apodization

%spectrum=abs(FourierDir(t,interferogram'.*apod,f));
spectrum=abs(FourierDir(t(pos_min:pos_max),interferogram(pos_min:pos_max).*apod,f));

spectrum_y_fixed=0; %condition 'do not fix spectrum y-axis'
set(handles.FixYSpectrum,'Value',spectrum_y_fixed);

interf_y_fixed=0; %condition 'do not fix interferogram y-axis'
set(handles.FixYInterferogram,'Value',interf_y_fixed);

plot(handles.Plot_spectr,f,spectrum,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_spectr,'XLim',[min(f) max(f)]);

if spectrum_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_spectr,'YLim',spectrum_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    spectrum_y_lim=get(handles.Plot_spectr,'YLim'); %take the new limits
end

plot(handles.Plot_Interf,t,interferogram,'Color',[0 0.7 0],'linewidth',0.5);
hold off;
set(handles.Plot_Interf,'XLim',[t_min t_max]);

if interf_y_fixed==1 %if the y axis is fixed
    set(handles.Plot_Interf,'YLim',interf_y_lim); %keep the y axis limits
else %if the y axis is not fixed
    interf_y_lim=get(handles.Plot_Interf,'YLim'); %take the new limits
end

% --- Executes during object creation, after setting all properties.
function t_maxValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_maxValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
