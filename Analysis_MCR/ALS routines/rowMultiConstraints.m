function varargout = rowMultiConstraints(varargin)
% ROWMULTICONSTRAINTS MATLAB code for rowMultiConstraints.fig
%      ROWMULTICONSTRAINTS, by itself, creates a new ROWMULTICONSTRAINTS or raises the existing
%      singleton*.
%
%      H = ROWMULTICONSTRAINTS returns the handle to a new ROWMULTICONSTRAINTS or the handle to
%      the existing singleton*.
%
%      ROWMULTICONSTRAINTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROWMULTICONSTRAINTS.M with the given input arguments.
%
%      ROWMULTICONSTRAINTS('Property','Value',...) creates a new ROWMULTICONSTRAINTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rowMultiConstraints_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rowMultiConstraints_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rowMultiConstraints

% Last Modified by GUIDE v2.5 04-Dec-2012 19:03:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @rowMultiConstraints_OpeningFcn, ...
    'gui_OutputFcn',  @rowMultiConstraints_OutputFcn, ...
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


% --- Executes just before rowMultiConstraints is made visible.
function rowMultiConstraints_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rowMultiConstraints (see VARARGIN)

% Choose default command line output for rowMultiConstraints
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
numbExp=evalin('base','mcr_als.alsOptions.nexp');
if numbExp==1
    set(handles.push_multiway,'enable','off');
else
    set(handles.push_multiway,'enable','on');
end

evalin('base','mcr_als.alsOptions.multi.ccons=0;'); % equal constraint to all Cs
evalin('base','mcr_als.alsOptions.multi.curr_cmat=1;'); % current matrix

% default isp
matc=evalin('base','mcr_als.alsOptions.multi.matc');
nsp=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
isp=ones(matc,nsp);
assignin('base','isp',isp);
evalin('base','mcr_als.alsOptions.multi.isp=isp;');
evalin('base','clear isp');

% total Nr of row submetrics
set(handles.text_num,'String',num2str(matc));

% UIWAIT makes rowMultiConstraints wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = rowMultiConstraints_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% *********************************************************************
% Actual Matrix

% --- Executes on button press in check_sameC.
function check_sameC_Callback(hObject, eventdata, handles)
% hObject    handle to check_sameC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_sameC

ccons=get(hObject,'Value');

if ccons==0
    set(handles.text_MatrixNr,'enable','on');
    matrixNR=evalin('base','mcr_als.alsOptions.multi.matc');
    j=1;
    for i=1:matrixNR,
        lsb=[i];
        lsv(j)={lsb};
        j=j+1;
    end;
    set(handles.popup_matrixNr,'enable','on','string',lsv,'value',1);
    
    %reset all values
    set(handles.check_nnC,'enable','on','value',0);
    set(handles.text_nnCi,'enable','off');
    set(handles.popup_nnCi,'enable','off','value',1);
    set(handles.text_nnCsp,'enable','off');
    set(handles.popup_nnCsp,'enable','off','value',1);
    set(handles.text_nnCp,'enable','off');
    set(handles.edit_nnCp,'enable','off','string',' ');
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''nonegC'');');
    
    set(handles.check_unimod,'enable','on','value',0);
    set(handles.text_Ui,'enable','off');
    set(handles.popup_Ui,'enable','off','value',1);
    set(handles.text_Up,'enable','off');
    set(handles.popup_Up,'enable','off','value',1);
    set(handles.text_Utol,'enable','off');
    set(handles.edit_Utol,'enable','off','string',' ');
    set(handles.text_Uv,'enable','off');
    set(handles.edit_Uv,'enable','off','string',' ');
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''unimodC'');');
    
    set(handles.check_closC,'enable','on','value',0);
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
    
    set(handles.check_eqC,'enable','on','value',0);
    set(handles.text_eqCs,'enable','off');
    set(handles.popup_eqCs,'enable','off','value',1);
    set(handles.text_eqCc,'enable','off');
    set(handles.popup_eqCc,'enable','off','value',1);
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''cselcC'');');
    
    evalin('base','mcr_als.alsOptions.nonegC.noneg=0;');
    evalin('base','mcr_als.alsOptions.unimodC.unimodal=0;');
    evalin('base','mcr_als.alsOptions.closure.closure=0;');
    evalin('base','mcr_als.alsOptions.cselcC.cselcon=0;');
    
else
    set(handles.text_MatrixNr,'enable','off');
    set(handles.popup_matrixNr,'enable','off','Value',1,'String','Same constraints');
    
    %reset all values
    set(handles.check_nnC,'enable','on','value',0);
    set(handles.text_nnCi,'enable','off');
    set(handles.popup_nnCi,'enable','off','value',1);
    set(handles.text_nnCsp,'enable','off');
    set(handles.popup_nnCsp,'enable','off','value',1);
    set(handles.text_nnCp,'enable','off');
    set(handles.edit_nnCp,'enable','off','string',' ');
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''nonegC'');');
    
    set(handles.check_unimod,'enable','on','value',0);
    set(handles.text_Ui,'enable','off');
    set(handles.popup_Ui,'enable','off','value',1);
    set(handles.text_Up,'enable','off');
    set(handles.popup_Up,'enable','off','value',1);
    set(handles.text_Utol,'enable','off');
    set(handles.edit_Utol,'enable','off','string',' ');
    set(handles.text_Uv,'enable','off');
    set(handles.edit_Uv,'enable','off','string',' ');
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''unimodC'');');
    
    set(handles.check_closC,'enable','on','value',0);
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
    
    set(handles.check_eqC,'enable','on','value',0);
    set(handles.text_eqCs,'enable','off');
    set(handles.popup_eqCs,'enable','off','value',1);
    set(handles.text_eqCc,'enable','off');
    set(handles.popup_eqCc,'enable','off','value',1);
    evalin('base','mcr_als.alsOptions=rmfield(mcr_als.alsOptions,''cselcC'');');
    
    evalin('base','mcr_als.alsOptions.nonegC.noneg=0;');
    evalin('base','mcr_als.alsOptions.unimodC.unimodal=0;');
    evalin('base','mcr_als.alsOptions.closure.closure=0;');
    evalin('base','mcr_als.alsOptions.cselcC.cselcon=0;');
    
end

assignin('base','ccons',ccons);
evalin('base','mcr_als.alsOptions.multi.ccons=ccons;');
evalin('base','clear ccons');


% --- Executes during object creation, after setting all properties.
function popup_matrixNr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_matrixNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

matrixNR=evalin('base','mcr_als.alsOptions.multi.matc');
j=1;
for i=1:matrixNR,
    lsb=[i];
    lsv(j)={lsb};
    j=j+1;
end;
set(hObject,'string',lsv)


% --- Executes on selection change in popup_matrixNr.
function popup_matrixNr_Callback(hObject, eventdata, handles)
% hObject    handle to popup_matrixNr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_matrixNr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_matrixNr

nmat = get(hObject,'value');
former_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');

curr_cmat=nmat;
assignin('base','curr_cmat',curr_cmat);
evalin('base','mcr_als.alsOptions.multi.curr_cmat=curr_cmat;');
evalin('base','clear curr_cmat');

% start -> values that could be modified

if get(handles.popup_nnCsp,'value') ~=1 & (former_cmat~=curr_cmat);
    treated=evalin('base','mcr_als.alsOptions.nonegC.treated');
    cneg=evalin('base','mcr_als.alsOptions.nonegC.cneg');
    if treated(1,3)==0
        % forzed to zero
        if treated(curr_cmat,1)==0
            set(handles.popup_nnCsp,'value',1);
            set(handles.edit_nnCp,'string',' ','enable','off');
        else
            valor=treated(curr_cmat,2)+2;
            nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
            nsign=nsign+2;
            set(handles.popup_nnCsp,'value',valor);
            if valor==(0+2)
                set(handles.edit_nnCp,'string',' ','enable','off');
            elseif valor==(nsign)
                set(handles.edit_nnCp,'string',' ','enable','off');
            else
                set(handles.edit_nnCp,'string',num2str(cneg(curr_cmat,:)),'enable','off');
            end
            
        end
    elseif treated(1,3)==1 | treated(1,3)==2
        % nnls and fnnls
        if treated(curr_cmat,1)==0
            set(handles.popup_nnCsp,'value',1);
            set(handles.edit_nnCp,'string',' ','enable','off');
        else
            valor=treated(curr_cmat,2)+2;
            set(handles.popup_nnCsp,'value',valor);
            set(handles.edit_nnCp,'string',' ','enable','off');
        end
    end
end

if get(handles.popup_Up,'value') ~=1 & (former_cmat~=curr_cmat);
    treated=evalin('base','mcr_als.alsOptions.unimodC.treated');
    spmod=evalin('base','mcr_als.alsOptions.unimodC.spmod');
    if treated(curr_cmat,1)==0
        set(handles.popup_Up,'value',1);
        set(handles.edit_Uv,'string',' ','enable','off');
    else
        valor=treated(curr_cmat,2)+2;
        nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
        nsign=nsign+2;
        set(handles.popup_Up,'value',valor);
        if valor==(0+2)
            set(handles.edit_Uv,'string',' ','enable','off');
        elseif valor==(nsign)
            set(handles.edit_Uv,'string',' ','enable','off');
        else
            set(handles.edit_Uv,'string',num2str(spmod(curr_cmat,:)),'enable','off');
        end
        
    end
end

if get(handles.popup_Cnr,'value') ~=1 & (former_cmat~=curr_cmat);
    treated=evalin('base','mcr_als.alsOptions.closure.treated');
    sclos1=evalin('base','mcr_als.alsOptions.closure.sclos1');
    sclos2=evalin('base','mcr_als.alsOptions.closure.sclos2');
    
    if treated(curr_cmat,1)==0
        set(handles.popup_Cnr,'value',1);
        set(handles.edit_C1,'string',' ','enable','off');
        set(handles.edit_C2,'string',' ','enable','off');
        set(handles.edit_C1v,'string',' ','enable','off');
        set(handles.edit_C2v,'string',' ','enable','off');
        set(handles.popup_C1c,'value',1,'enable','off');
        set(handles.popup_C2c,'value',1,'enable','off');
        set(handles.edit_C1sp,'string',' ','enable','off');
        set(handles.edit_C2sp,'string',' ','enable','off');
        set(handles.check_C1sp,'value',0,'enable','off');
        set(handles.check_C2sp,'value',0,'enable','off');
        set(handles.check_Cvar,'value',0,'enable','off');
    else
        % number of closures
        clos_number=treated(curr_cmat,2);
        if clos_number==0
            set(handles.popup_Cnr,'value',clos_number+2);
            set(handles.edit_C1,'string',' ','enable','off');
            set(handles.edit_C2,'string',' ','enable','off');
            set(handles.edit_C1v,'string',' ','enable','off');
            set(handles.edit_C2v,'string',' ','enable','off');
            set(handles.popup_C1c,'value',1,'enable','off');
            set(handles.popup_C2c,'value',1,'enable','off');
            set(handles.edit_C1sp,'string',' ','enable','off');
            set(handles.edit_C2sp,'string',' ','enable','off');
            set(handles.check_C1sp,'value',0,'enable','off');
            set(handles.check_C2sp,'value',0,'enable','off');
            set(handles.check_Cvar,'value',0,'enable','off');
        elseif clos_number==1
            set(handles.popup_Cnr,'value',clos_number+2);
            set(handles.edit_C1,'string',num2str(treated(curr_cmat,3)),'enable','on');
            set(handles.edit_C2,'string',' ','enable','off');
            set(handles.edit_C1v,'string',' ','enable','off');
            set(handles.edit_C2v,'string',' ','enable','off');
            set(handles.popup_C1c,'value',treated(curr_cmat,4)+1,'enable','on');
            set(handles.popup_C2c,'value',1,'enable','off');
            set(handles.edit_C1sp,'string',num2str(sclos1(curr_cmat,:)),'enable','on');
            set(handles.edit_C2sp,'string',' ','enable','off');
            
            nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
            clos_sp1=treated(curr_cmat,5);
            if nsign==clos_sp1
                set(handles.check_C1sp,'value',1,'enable','on');
            else
                set(handles.check_C1sp,'value',0,'enable','on');
            end
            set(handles.check_C2sp,'value',0,'enable','off');
            
            set(handles.check_Cvar,'value',0,'enable','off');
            
        elseif clos_number==2
            set(handles.popup_Cnr,'value',clos_number+2);
            set(handles.edit_C1,'string',num2str(treated(curr_cmat,3)),'enable','on');
            set(handles.edit_C2,'string',num2str(treated(curr_cmat,6)),'enable','on');
            set(handles.edit_C1v,'string',' ','enable','off');
            set(handles.edit_C2v,'string',' ','enable','off');
            set(handles.popup_C1c,'value',treated(curr_cmat,4)+1,'enable','on');
            set(handles.popup_C2c,'value',treated(curr_cmat,7)+1,'enable','on');
            set(handles.edit_C1sp,'string',num2str(sclos1(curr_cmat,:)),'enable','on');
            set(handles.edit_C2sp,'string',num2str(sclos2(curr_cmat,:)),'enable','on');
            
            nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
            clos_sp1=treated(curr_cmat,5);
            clos_sp2=treated(curr_cmat,8);
            if nsign==clos_sp1
                set(handles.check_C1sp,'value',1,'enable','on');
            else
                set(handles.check_C1sp,'value',0,'enable','on');
            end
            if nsign==clos_sp2
                set(handles.check_C2sp,'value',1,'enable','on');
            else
                set(handles.check_C2sp,'value',0,'enable','on');
            end
            
            set(handles.check_Cvar,'value',0,'enable','off');
            
        end
        
    end
end


% end -> values that could be modified

assignin('base','curr_cmat',curr_cmat);
evalin('base','mcr_als.alsOptions.multi.curr_cmat=curr_cmat;');
evalin('base','clear curr_cmat');

% *********************************************************************
% ISP

% --- Executes during object creation, after setting all properties.
function popup_isp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_isp (see GCBO)
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
lsv(1)={'select...'};
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

% --- Executes on selection change in popup_isp.
function popup_isp_Callback(hObject, eventdata, handles)
% hObject    handle to popup_isp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_isp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_isp

popmenu3=get(handles.popup_isp,'String');
pm3=get(handles.popup_isp,'Value');

if pm3==1
    
else
    ispmat_char=char([popmenu3(pm3)]);
    ispmat=evalin('base',ispmat_char);
    
    assignin('base','isp',ispmat);
    assignin('base','ispmat_char',ispmat_char);
    
    evalin('base','mcr_als.alsOptions.multi.isp=isp;');
    evalin('base','mcr_als.alsOptions.multi.ispmat_char=ispmat_char;');
    
    evalin('base','clear isp ispmat_char');
    
end


% --- Executes on button press in check_ISPdefault.
function check_ISPdefault_Callback(hObject, eventdata, handles)
% hObject    handle to check_ISPdefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_ISPdefault

ispcon=get(hObject,'Value');

if ispcon==1;
    matc=evalin('base','mcr_als.alsOptions.multi.matc');
    nsp=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    isp=ones(matc,nsp);
    
    assignin('base','isp',isp);
    evalin('base','mcr_als.alsOptions.multi.isp=isp;');
    evalin('base','clear isp');
    
    
    set(handles.popup_isp,'enable','off','value',1);
    set(handles.text_popISP,'enable','off');
    
else
    
    set(handles.popup_isp,'enable','on','value',1);
    set(handles.text_popISP,'enable','on');
    
end;



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
    evalin('base','mcr_als.alsOptions.nonegC.cneg=[];');
    
    evalin('base','mcr_als.alsOptions.nonegC.treated=zeros(mcr_als.alsOptions.multi.matc,3);');
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
    
    dim=evalin('base','min(size(mcr_als.alsOptions.iniesta));');
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
    dim=evalin('base','min(size(mcr_als.alsOptions.iniesta));');
    llistannc(1)={'select...'};
    llistnnc=[0 dim];
    llistannc(2)={llistnnc};
    
    set(handles.text_nnCsp,'enable','on');
    set(handles.popup_nnCsp,'enable','on','value',1)
    set(handles.popup_nnCsp,'string',llistannc)
    
    set(handles.text_nnCp,'enable','off');
    set(handles.edit_nnCp,'enable','off','string','');
    
end

evalin('base','mcr_als.alsOptions.nonegC.treated(:,3)=ialg;');

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
matc=evalin('base','mcr_als.alsOptions.multi.matc');
nexp=evalin('base','mcr_als.alsOptions.nexp');
ialg=evalin('base','mcr_als.alsOptions.nonegC.ialg');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
cneg=evalin('base','mcr_als.alsOptions.nonegC.cneg');
assignin('base','matc',matc);
assignin('base','dim',dim);
assignin('base','curr_cmat',curr_cmat);

if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    if ialg == 0
        if ncneg == dim;
            cneg = ones(matc,dim);
            set(handles.text_nnCp,'enable','off');
            set(handles.edit_nnCp,'enable','off','string','');
        elseif  ncneg == -1
            cneg=zeros(matc,dim);
            set(handles.text_nnCp,'enable','off');
            set(handles.edit_nnCp,'enable','off','string','');
        elseif  ncneg == 0
            cneg = zeros(matc,dim);
            set(handles.text_nnCp,'enable','off');
            set(handles.edit_nnCp,'enable','off','string','');
        else
            cneg=[];
            set(handles.text_nnCp,'enable','on');
            set(handles.edit_nnCp,'enable','on');
        end
    else
        %         cneg = ones(nexp,dim);
        if ncneg==0;
            cneg = zeros(matc,dim);
        elseif ncneg==1;
            cneg = ones(matc,dim);
        end
        
    end
elseif evalin('base','mcr_als.alsOptions.multi.ccons')==0
    if ialg == 0
        if ncneg == dim;
            cnegr=ones(1,dim);
            cneg(curr_cmat,:)=cnegr;
            assignin('base','cneg',cneg);
            evalin('base','mcr_als.alsOptions.nonegC.cneg = cneg;');
            set(handles.text_nnCp,'enable','off');
            set(handles.edit_nnCp,'enable','off','string','');
        elseif  ncneg == -1
            cnegr=zeros(1,dim);
            cneg(curr_cmat,:)=cnegr;
            assignin('base','cneg',cneg);
            evalin('base','mcr_als.alsOptions.nonegC.cneg = cneg;');
            set(handles.text_nnCp,'enable','off');
            set(handles.edit_nnCp,'enable','off','string','');
        elseif  ncneg == 0
            cnegr=zeros(1,dim);
            cneg(curr_cmat,:)=cnegr;
            assignin('base','cneg',cneg);
            evalin('base','mcr_als.alsOptions.nonegC.cneg = cneg;');
            set(handles.text_nnCp,'enable','off');
            set(handles.edit_nnCp,'enable','off','string','');
        else
            set(handles.text_nnCp,'enable','on');
            set(handles.edit_nnCp,'enable','on');
        end
    else
        if ncneg==0;
            form_cneg = zeros(1,dim);
            cneg(curr_cmat,:)=form_cneg;
        elseif ncneg==1;
            form_cneg= ones(1,dim);
            cneg(curr_cmat,:)=form_cneg;
        end
    end
    
    evalin('base','mcr_als.alsOptions.nonegC.treated(curr_cmat,1)=1;');
    evalin('base','mcr_als.alsOptions.nonegC.treated(curr_cmat,2)=ncneg;');
    
end


assignin('base','cneg',cneg);
evalin('base','mcr_als.alsOptions.nonegC.cneg=cneg;');
evalin('base','clear cneg ncneg curr_cmat matc dim');


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

curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
cneg=evalin('base','mcr_als.alsOptions.nonegC.cneg');
assignin('base','curr_cmat',curr_cmat);
cneg(curr_cmat,:)=cnegr;

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
    
    evalin('base','mcr_als.alsOptions.unimodC.spmod=[];');
    % default value of Utol
    rmod=1.1;
    assignin('base','rmod',rmod);
    evalin('base','mcr_als.alsOptions.unimodC.rmod=rmod;');
    evalin('base','clear rmod');
    
    evalin('base','mcr_als.alsOptions.unimodC.treated=zeros(mcr_als.alsOptions.multi.matc,2);');
    
    
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

assignin('base','nmod',nmod);
evalin('base','mcr_als.alsOptions.unimodC.nmod=nmod;');

dim=evalin('base','min(size(mcr_als.alsOptions.iniesta));');
nexp=evalin('base','mcr_als.alsOptions.nexp;');
matc=evalin('base','mcr_als.alsOptions.multi.matc;');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat;');
spmod=evalin('base','mcr_als.alsOptions.unimodC.spmod;');
assignin('base','matc',matc);
assignin('base','dim',dim);
assignin('base','curr_cmat',curr_cmat);


if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    
    if nmod == dim
        spmod=ones(matc,dim);
        set(handles.text_Uv,'enable','off');
        set(handles.edit_Uv,'enable','off','string','');
    elseif  nmod == 0
        spmod=zeros(matc,dim);
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
    
    
elseif evalin('base','mcr_als.alsOptions.multi.ccons')==0
    
    if nmod == dim
        form_spmod=ones(1,dim);
        spmod(curr_cmat,:)=form_spmod;
        set(handles.text_Uv,'enable','off');
        set(handles.edit_Uv,'enable','off','string','');
    elseif  nmod == 0
        form_spmod=zeros(1,dim);
        spmod(curr_cmat,:)=form_spmod;
        set(handles.text_Uv,'enable','off');
        set(handles.edit_Uv,'enable','off','string','');
    elseif  nmod == -1
        evalin('base','mcr_als.alsOptions.unimodC.spmod=[]');
        set(handles.text_Uv,'enable','off');
        set(handles.edit_Uv,'enable','off','string','');
    else
        evalin('base','mcr_als.alsOptions.unimodC.spmod=[]');
        set(handles.text_Uv,'enable','on');
        set(handles.edit_Uv,'enable','on');
    end
    
    evalin('base','mcr_als.alsOptions.unimodC.treated(curr_cmat,1)=1;');
    evalin('base','mcr_als.alsOptions.unimodC.treated(curr_cmat,2)=nmod;');
    
    
end

assignin('base','spmod',spmod);
evalin('base','mcr_als.alsOptions.unimodC.spmod=spmod;');
evalin('base','clear spmod curr_cmat dim matc nmod');



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

form_spmod= str2num(get(hObject,'String'));
matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
spmod=evalin('base','mcr_als.alsOptions.unimodC.spmod');

if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    for i=1:matc;
        spmod(i,:)=form_spmod;
    end
elseif evalin('base','mcr_als.alsOptions.multi.ccons')==0
    spmod(curr_cmat,:)=form_spmod;
end

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
    
    dim=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    assignin('base','dim',dim);
    matc=evalin('base','mcr_als.alsOptions.multi.matc');
    assignin('base','matc',matc);
    evalin('base','mcr_als.alsOptions.closure.tclos1(1:matc)=zeros(1,matc);');
    evalin('base','mcr_als.alsOptions.closure.tclos2(1:matc)=zeros(1,matc);');
    evalin('base','mcr_als.alsOptions.closure.sclos1(1:matc,1:dim)=zeros(matc,dim);');
    evalin('base','mcr_als.alsOptions.closure.sclos2(1:matc,1:dim)=zeros(matc,dim);');
    evalin('base','mcr_als.alsOptions.closure.iclos(1:matc)=zeros(1,matc);');
    evalin('base','mcr_als.alsOptions.closure.iclos1(1:matc)=zeros(1,matc);');
    evalin('base','mcr_als.alsOptions.closure.iclos2(1:matc)=zeros(1,matc);');
    evalin('base','mcr_als.alsOptions.closure.vclos1(1:matc)=zeros(1,matc);');
    evalin('base','mcr_als.alsOptions.closure.vclos2(1:matc)=zeros(1,matc);');
    
    evalin('base','mcr_als.alsOptions.closure.treated=zeros(matc,8);'); % treated/nr closures/valor 1/condition1 /nr species 1/valor 2/condition 2/ nr species 2
    evalin('base','clear dc dim matc');
    
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
llista(2)={'0'};
llista(3)={'1'};
llista(4)={'2'};
set(hObject,'string',llista)

% --- Executes on selection change in popup_Cnr.
function popup_Cnr_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Cnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Cnr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Cnr

form_iclos=get(hObject,'Value')-2;

matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
assignin('base','matc',matc);
assignin('base','curr_cmat',curr_cmat);

if form_iclos == 1;
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
    
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,1)=1;');
    
    
elseif  form_iclos==2
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
    
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,1)=1;');
    
elseif form_iclos == 0;
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
    
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,1)=1;');
    
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
    
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,1)=0;');
end

% matr=evalin('base','mcr_als.alsOptions.multi.matr');
matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
% curr_smat=evalin('base','mcr_als.alsOptions.multi.curr_smat');
iclos=evalin('base','mcr_als.alsOptions.closure.iclos');

if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        iclos(1,[1:matc])=form_iclos;
    end
    
else
    
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        iclos(1,curr_cmat)=form_iclos;
    end
    
end

assignin('base','iclos',iclos);
evalin('base','mcr_als.alsOptions.closure.iclos=iclos;');
assignin('base','form_iclos',form_iclos);
evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,2)=form_iclos;');
evalin('base','clear iclos curr_cmat matc');


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

form_tclos1= str2num(get(hObject,'String'));

matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
assignin('base','matc',matc);
assignin('base','curr_cmat',curr_cmat);

tclos1=evalin('base','mcr_als.alsOptions.closure.tclos1');

if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        tclos1(1,[1:matc])=form_tclos1;
    end
    
else
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        tclos1(1,curr_cmat)=form_tclos1;
    end
end

assignin('base','tclos1',tclos1);
evalin('base','mcr_als.alsOptions.closure.tclos1=tclos1;');

assignin('base','form_tclos1',form_tclos1);
evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,3)=form_tclos1;');
evalin('base','clear tclos1 matc curr_cmat form_tclos1');


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

form_iclos1=get(hObject,'Value')-1;

matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
assignin('base','matc',matc);
assignin('base','curr_cmat',curr_cmat);

iclos1=evalin('base','mcr_als.alsOptions.closure.iclos1');

if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        iclos1(1,[1:matc])=form_iclos1;
    end
    
else
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        iclos1(1,curr_cmat)=form_iclos1;
        
    end
end

assignin('base','iclos1',iclos1);
evalin('base','mcr_als.alsOptions.closure.iclos1=iclos1;');

assignin('base','form_iclos1',form_iclos1);
evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,4)=form_iclos1;');
evalin('base','clear iclos1 matc curr_cmat form_iclos1');


% --- Executes on button press in check_C1sp.
function check_C1sp_Callback(hObject, eventdata, handles)
% hObject    handle to check_C1sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_C1sp

check1clos=get(hObject,'Value');

matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
assignin('base','matc',matc);
assignin('base','curr_cmat',curr_cmat);

if check1clos==0
    set(handles.text_C1sp,'enable','on');
    set(handles.edit_C1sp,'enable','on','string',' ');
    nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    ceros=zeros(1,nsign);
    sclos1=evalin('base','mcr_als.alsOptions.closure.sclos1');
    
    if evalin('base','mcr_als.alsOptions.multi.ccons')==1
        if evalin('base','mcr_als.alsOptions.closure.dc')==1
            for i=1:matc;
                sclos1(i,:)=ceros;
            end
        end
    else
        if evalin('base','mcr_als.alsOptions.closure.dc')==1
            sclos1(curr_cmat,:)=ceros;
        end
    end
    
    assignin('base','sclos1',sclos1);
    evalin('base','mcr_als.alsOptions.closure.sclos1=sclos1;');
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,5)=0;');
    evalin('base','clear sclos1 matc curr_cmat');
    
    
elseif check1clos==1
    nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    unos=ones(1,nsign);
    set(handles.text_C1sp,'enable','on');
    set(handles.edit_C1sp,'enable','inactive','string',num2str(unos));
    
    sclos1=evalin('base','mcr_als.alsOptions.closure.sclos1');
    
    if evalin('base','mcr_als.alsOptions.multi.ccons')==1
        if evalin('base','mcr_als.alsOptions.closure.dc')==1
            for i=1:matc;
                sclos1(i,:)=unos;
            end
        end
    else
        if evalin('base','mcr_als.alsOptions.closure.dc')==1
            sclos1(curr_cmat,:)=unos;
        end
    end
    
    assignin('base','sclos1',sclos1);
    evalin('base','mcr_als.alsOptions.closure.sclos1=sclos1;');
    suma_sp=sum(sclos1(curr_cmat,:));
    assignin('base','suma_sp',suma_sp);
    
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,5)=suma_sp;');
    evalin('base','clear sclos1 matc curr_cmat suma_sp');
    
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


form_sclos1= str2num(get(hObject,'String'));

matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
assignin('base','matc',matc);
assignin('base','curr_cmat',curr_cmat);
sclos1=evalin('base','mcr_als.alsOptions.closure.sclos1');

if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        for i=1:matc;
            sclos1(i,:)=form_sclos1;
        end
    end
else
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        sclos1(curr_cmat,:)=form_sclos1;
    end
end

assignin('base','sclos1',sclos1);
evalin('base','mcr_als.alsOptions.closure.sclos1=sclos1;');
suma_sp=sum(form_sclos1);
assignin('base','suma_sp',suma_sp);
evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,5)=suma_sp;');
evalin('base','clear sclos1 matc curr_cmat suma_sp');


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

form_tclos2= str2num(get(hObject,'String'));

matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
assignin('base','matc',matc);
assignin('base','curr_cmat',curr_cmat);

tclos2=evalin('base','mcr_als.alsOptions.closure.tclos2');

if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        tclos2(1,[1:matc])=form_tclos2;
    end
    
else
    
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        tclos2(1,curr_cmat)=form_tclos2;
    end
    
end

assignin('base','tclos2',tclos2);
evalin('base','mcr_als.alsOptions.closure.tclos2=tclos2;');

assignin('base','form_tclos2',form_tclos2);
evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,6)=form_tclos2;');

evalin('base','clear tclos2 matc curr_cmat form_tclos2');



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


form_iclos2=get(hObject,'Value')-1;

matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
assignin('base','matc',matc);
assignin('base','curr_cmat',curr_cmat);

iclos2=evalin('base','mcr_als.alsOptions.closure.iclos1');

if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        iclos2(1,[1:matc])=form_iclos2;
    end
    
    
else
    
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        iclos2(1,curr_cmat)=form_iclos2;
    end
    
end

assignin('base','iclos2',iclos2);
evalin('base','mcr_als.alsOptions.closure.iclos2=iclos2;');

assignin('base','form_iclos2',form_iclos2);
evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,7)=form_iclos2;');
evalin('base','clear iclos2 matc curr_cmat form_iclos2');


% --- Executes on button press in check_C2sp.
function check_C2sp_Callback(hObject, eventdata, handles)
% hObject    handle to check_C2sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_C2sp


check2clos=get(hObject,'Value');

matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
assignin('base','matc',matc);
assignin('base','curr_cmat',curr_cmat);

if check2clos==0
    nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    ceros=zeros(1,nsign);
    set(handles.text_C2sp,'enable','on');
    set(handles.edit_C2sp,'enable','on','string',' ');
    
    sclos2=evalin('base','mcr_als.alsOptions.closure.sclos2');
    
    if evalin('base','mcr_als.alsOptions.multi.ccons')==1
        if evalin('base','mcr_als.alsOptions.closure.dc')==1
            for i=1:matc;
                sclos2(i,:)=ceros;
            end
            
        end
    else
        if evalin('base','mcr_als.alsOptions.closure.dc')==1
            sclos2(curr_cmat,:)=ceros;
        end
    end
    
    assignin('base','sclos2',sclos2);
    evalin('base','mcr_als.alsOptions.closure.sclos1=sclos2;');
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,8)=0')
    evalin('base','clear sclos2 matc curr_cmat');
    
elseif check2clos==1
    nsign=evalin('base','min(size(mcr_als.alsOptions.iniesta))');
    unos=ones(1,nsign);
    set(handles.text_C2sp,'enable','on');
    set(handles.edit_C2sp,'enable','inactive','string',num2str(unos));
    
    sclos2=evalin('base','mcr_als.alsOptions.closure.sclos2');
    
    if evalin('base','mcr_als.alsOptions.multi.ccons')==1
        if evalin('base','mcr_als.alsOptions.closure.dc')==1
            for i=1:matc;
                sclos2(i,:)=unos;
            end
        end
    else
        if evalin('base','mcr_als.alsOptions.closure.dc')==1
            sclos2(curr_cmat,:)=unos;
        end
    end
    
    assignin('base','sclos2',sclos2);
    evalin('base','mcr_als.alsOptions.closure.sclos2=sclos2;');
    
    suma_sp=sum(sclos2(curr_cmat,:));
    assignin('base','suma_sp',suma_sp);
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,8)=suma_sp;');
    evalin('base','clear sclos2 matc curr_cmat suma_sp');
    
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

form_sclos2= str2num(get(hObject,'String'));

matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
assignin('base','matc',matc);
assignin('base','curr_cmat',curr_cmat);

sclos2=evalin('base','mcr_als.alsOptions.closure.sclos2');

if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        for i=1:matc;
            sclos2(i,:)=form_sclos2;
        end
    end
    
else
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        sclos2(curr_cmat,:)=form_sclos2;
    end
end

assignin('base','sclos2',sclos2);
evalin('base','mcr_als.alsOptions.closure.sclos2=sclos2;');

suma_sp=sum(form_sclos2);
assignin('base','suma_sp',suma_sp);
evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,8)=suma_sp;');

evalin('base','clear sclos2 matc curr_cmat form_sclos2 suma_sp');


% vclos
% oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


% --- Executes on button press in check_Cvar.
function check_Cvar_Callback(hObject, eventdata, handles)
% hObject    handle to check_Cvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_Cvar


vclos=get(hObject,'Value');

matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
assignin('base','matc',matc);
assignin('base','curr_cmat',curr_cmat);


if vclos==1;
    iclos=evalin('base','mcr_als.alsOptions.closure.iclos');
    
    if evalin('base','mcr_als.alsOptions.closure.dc')==1;
        activa=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
    elseif evalin('base','mcr_als.alsOptions.closure.dc')==2;
        activa=evalin('base','mcr_als.alsOptions.multi.curr_smat');
    end
    
    if iclos(activa)==1
        
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,3)=Inf');
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,4)=0;');
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,5)=0;');
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,6)=0;');
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,7)=0;');
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,8)=0;');
        
        set(handles.text_C1v,'enable','on');
        set(handles.edit_C1v,'enable','on','string',' ');
        set(handles.text_C2v,'enable','off');
        set(handles.edit_C2v,'enable','off','string',' ');
        set(handles.text_C1,'enable','off');
        set(handles.edit_C1,'enable','off','string',' ');
        set(handles.text_C2,'enable','off');
        set(handles.edit_C2,'enable','off','string',' ');
        set(handles.text_C1c,'enable','off');
        set(handles.popup_C1c,'enable','off','value',1);
        set(handles.text_C2c,'enable','off');
        set(handles.popup_C2c,'enable','off','value',1);
        set(handles.text_C1sp,'enable','on');
        set(handles.edit_C1sp,'enable','on','value',1);
        set(handles.text_C2sp,'enable','off');
        set(handles.edit_C2sp,'enable','off','value',1);
        
        iclos1=1; % equal
        assignin('base','iclos1',iclos1);
        evalin('base','mcr_als.alsOptions.closure.iclos1=iclos1;');
        evalin('base','clear iclos1');
        
        vc=1;  % equal
        assignin('base','vc',vc);
        evalin('base','mcr_als.alsOptions.closure.vc=vc;');
        evalin('base','clear vc');
        
    elseif iclos(activa)==2
        
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,3)=Inf');
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,4)=0;');
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,5)=0;');
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,6)=Inf;');
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,7)=0;');
        evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,8)=0;');
        
        set(handles.text_C1v,'enable','on');
        set(handles.edit_C1v,'enable','on','string',' ');
        set(handles.text_C2v,'enable','on');
        set(handles.edit_C2v,'enable','on','string',' ');
        set(handles.text_C1,'enable','off');
        set(handles.edit_C1,'enable','off','string',' ');
        set(handles.text_C2,'enable','off');
        set(handles.edit_C2,'enable','off','string',' ');
        set(handles.text_C1c,'enable','off');
        set(handles.popup_C1c,'enable','off','value',1);
        set(handles.text_C2c,'enable','off');
        set(handles.popup_C2c,'enable','off','value',1);
        set(handles.text_C1sp,'enable','on');
        set(handles.edit_C1sp,'enable','on','value',1);
        set(handles.text_C2sp,'enable','on');
        set(handles.edit_C2sp,'enable','on','value',1);
        
        iclos1=1; % equal
        assignin('base','iclos1',iclos1);
        evalin('base','mcr_als.alsOptions.closure.iclos1=iclos1;');
        evalin('base','clear iclos1');
        
        iclos2=1; % equal
        assignin('base','iclos2',iclos2);
        evalin('base','mcr_als.alsOptions.closure.iclos2=iclos2;');
        evalin('base','clear iclos2');
        
        vc=1; % equal
        assignin('base','vc',vc);
        evalin('base','mcr_als.alsOptions.closure.vc=vc;');
        evalin('base','clear vc');
        
    end
    
else
    
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,3)=0');
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,4)=0;');
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,5)=0;');
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,6)=0;');
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,7)=0;');
    evalin('base','mcr_als.alsOptions.closure.treated(curr_cmat,8)=0;');
    
    set(handles.text_C1v,'enable','off');
    set(handles.edit_C1v,'enable','off','string',' ');
    set(handles.text_C2v,'enable','off');
    set(handles.edit_C2v,'enable','off','string',' ');
    set(handles.text_Cnr,'enable','on');
    set(handles.popup_Cnr,'enable','on','value',1);
    set(handles.text_C1,'enable','off');
    set(handles.edit_C1,'enable','off','string',' ');
    set(handles.text_C2,'enable','off');
    set(handles.edit_C2,'enable','off','string',' ');
    set(handles.text_C1c,'enable','off');
    set(handles.popup_C1c,'enable','off','value',1);
    set(handles.text_C2c,'enable','off');
    set(handles.popup_C2c,'enable','off','value',1);
    set(handles.text_C1sp,'enable','off');
    set(handles.edit_C1sp,'enable','off','string',' ');
    set(handles.text_C2sp,'enable','off');
    set(handles.edit_C2sp,'enable','off','string',' ');
    
    vc=0; % equal
    assignin('base','vc',vc);
    evalin('base','mcr_als.alsOptions.closure.vc=vc;');
    evalin('base','clear vc');
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
form_vclos1=evalin('base',vclos1);

matc=evalin('base','mcr_als.alsOptions.multi.matc');

curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
vclos1=evalin('base','mcr_als.alsOptions.closure.vclos1');

if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        vclos1([1:matc],1)=form_vclos1;
    end
    
else
    
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        vclos1(curr_cmat,1)=form_vclos1;
    end
    
end

assignin('base','vclos1',vclos1);
evalin('base','mcr_als.alsOptions.closure.vclos1=vclos1;');
evalin('base','clear vclos1 matc curr_cmat');


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
form_vclos2=evalin('base',vclos2);

matc=evalin('base','mcr_als.alsOptions.multi.matc');
curr_cmat=evalin('base','mcr_als.alsOptions.multi.curr_cmat');
vclos2=evalin('base','mcr_als.alsOptions.closure.vclos2');

if evalin('base','mcr_als.alsOptions.multi.ccons')==1
    
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        vclos2([1:matc],1)=form_vclos2;
    end
    
else
    
    if evalin('base','mcr_als.alsOptions.closure.dc')==1
        vclos2(curr_cmat,1)=form_vclos2;
    end
    
end

assignin('base','vclos2',vclos2);
evalin('base','mcr_als.alsOptions.closure.vclos2=vclos2;');
evalin('base','clear vclos2 matc curr_cmat');


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

% Multiway
% ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

% --- Executes on button press in push_multiway.
function push_multiway_Callback(hObject, eventdata, handles)
% hObject    handle to push_multiway (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

appTril=evalin('base','mcr_als.alsOptions.trilin.appTril;');
multiwayC;
% if appTril==0
%     gMultiwayRow;
% else
%     gMultiway_reload;
% end


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

kineticHM_multi;


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
elseif evalin('base','mcr_als.alsOptions.multi.datamod')==3
    columnMultiConstraints;
end

% Reset
% *********************************************************************


% --- Executes on button press in push_reset.
function push_reset_Callback(hObject, eventdata, handles)
% hObject    handle to push_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close;
rowMultiConstraints;


% *********************************************************************
% *********************************************************************
% *********************************************************************


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
