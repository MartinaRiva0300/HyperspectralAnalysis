%Preparation of raw data for temporal hypercube generation- Martina
%15/10/2024

% --VERSIONE LINUX O WINDOWS-- 1
if isunix
    linux=1; %0 se si usa Windows, 1 (o altri valori) se si usa Linux
else
    linux=0;
end;

if ne(linux,0) %se si sta usando Linux
    fprintf('\n\n--LINUX VERSION of GetRawData, for data preparation--\n\n');
else %se non si sta usando linux
    fprintf('\n\n--WINDOWS VERSION of GetRawData, for data preparation--\n\n');
end

%**Benedetto**

% Grey colormap
M=linspace(0,1,512)';
MAP=[M M M];

%speed of light
c=299792458;

%**Benedetto**
if ne(linux,0) %se si sta usando Linux
    path(path,'/home/cristian/ProgrammiMatlab');
    path(path,'/home/cristian/ProgrammiMatlab/Wedges/Hyperspectral/RGBspectra');
    %dir0=('/home/cristian/Dropbox/Hyperspectral_image');

    % dir0=('/home/cristian/Work/Progetti/Brasile/Oct2017/BeloHorizonte/Microscope');

    dir0=('/home/cristian/Experiments/Wedges/Hyperspectral/Data/Microscope/');
    % dir0=('/home/cristian/Work/Experiments/Wedges/Hyperspectral/Data/Tests');
else %se non si sta usando linux
    path(path,'C:\Users\marti\OneDrive - Politecnico di Milano\MAGISTRALE\Thesis\Measurements');

    dir0=('C:\Users\marti\OneDrive - Politecnico di Milano\MAGISTRALE\Thesis\Measurements');
end


%Read file with motor positions

%save the .xls file in .txt with name "Delays" 
 dirTime=uigetdir(dir0, 'Select directory for raw data folder'); 
 txtTime= dir(fullfile(dirTime, '*.txt')); %txt file
 xlsTime=dir(fullfile(dirTime, '*.xls')); %xls file


 if not(isfolder([dirTime 'Infos']))
    mkdir('Infos');
 else 
     fprintf('\n Infos folder already exists in the chosen directory. Check your data!')
     return
 end;

newFolder=[dirTime, '\', 'Infos'];
mkdir(newFolder)
newFolderIn=[newFolder, '\', 'Old'];
mkdir(newFolderIn)
newFolderIn2=[newFolder, '\', 'Results'];
mkdir(newFolderIn2)


%moving files in Infos 
movefile([dirTime, '\', txtTime.name], [dirTime, '\', 'Infos', '\', 'Delays.txt'])
movefile([dirTime, '\', xlsTime.name], newFolderIn)

%rename bins
RenameBin

fprintf('\n Now your data are ready to be processed with HyperspectralAnalysis_Time')
