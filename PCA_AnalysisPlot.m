function PCA_AnalysisPlot(n_PCA,Mask_selection)
%This function plot the first n_PCA Principal Components (PCs) for the PCA
%of a selected Spectral Hypercube
%If Mask_selection is equal to 1 then the selection of a mask can be
%performed

MAT_HYP=0; %MAT_HYP=1 load .mat hypercube

if MAT_HYP==1
    dir0=('C:\Users\HARDi\Desktop\Tesi');
    [filename, pathname] = uigetfile('*.mat', 'Load Spectral Hypercube',dir0);
    file_tot=[pathname,filename];

    h=waitbar(1/2,'Loading the Spectra Hypercube');
    input_data=load(file_tot);
    close(h);
    Hyperspectrum_cube=abs(input_data.Hyperspectrum_cube);
    fr_real=input_data.fr_real;
    c=299792458; %speed of light [m/s]
    wavelength=c./fr_real.*1e-3;
    
    if isfield(input_data,'Intens')==1 %if Intens has been saved together with Hyperspectrum_cube **Benedetto** 12/07/2021
        Intens=input_data.Intens;
    else
        Intens=sum(abs(Hyperspectrum_cube),3);
    end
else
    [Hyperspectrum_cube,f,saturationMap,file_totCal,fr_real,Intens] = H5HypercubeRead;
    c=3*1e8;
    wavelength=c./fr_real.*1e-3;
    if isempty(Intens) %if Intens has not been saved together with Hyperspectrum_cube **Benedetto** 12/07/2021
        Intens=sum(abs(Hyperspectrum_cube),3);
    end
end

max_wavelength=max(wavelength);
min_wavelength=min(wavelength);

if Mask_selection==1
    %MASK SELECTION PART
    
    h_fig=figure;
    imagesc(Intens);
    M=linspace(0,1,512)';
    MAP=[M M M];
    colormap(MAP);
    title('Please select the AOI mask');

    [ThisAOI] = roipoly();

    for index=1:length(fr_real)
        Frame=Hyperspectrum_cube(:,:,index);
        Frame(ThisAOI==0)=NaN;
        Hyperspectrum_cube(:,:,index)=Frame;
    end
    
    Intens(ThisAOI==0)=NaN; %The Intens map is updated according to the mask
    
    close(h_fig);
end

%ROTATION SPECIFICATION PART
h_fig=figure;
imagesc(Intens);
M=linspace(0,1,512)';
MAP=[M M M];
colormap(MAP);

prompt={'Select rotation angle (+: clock, -: counterclock)'};
name='Select rotation angle';
numlines=1;
defaultanswer={num2str(0)};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer=inputdlg(prompt,name,numlines,defaultanswer,options);
angle=str2double(cell2mat(answer));
close(h_fig);
%END ROTATION SPECIFICATION PART

hw2=waitbar(1/2,'Hypercube PCA processing. Please wait...');

[COEFF,SCORE]=PCA_tool(Hyperspectrum_cube);

max_COEFF=max(max(COEFF(:,1:n_PCA)));
min_COEFF=min(min(COEFF(:,1:n_PCA)));


pos1 = [100 50 1300 700];
pos2 = [0.1 0.3 0.3 0.3];

for index=1:n_PCA
    figure;
    set(gcf,'Position',pos1)
    
    subplot('Position',pos2)
    subplot(1,5,[1,2]);
    plot(wavelength,COEFF(:,index),'Linewidth',2);
    xlabel('wavelength [nm]', 'FontSize', 16);
    ylabel('coefficient', 'FontSize', 16);
    t1=title(['PC ',num2str(index)]);
    t1.FontSize=16;
    ax = gca;
    ax.FontSize = 16; 
    axis([min_wavelength max_wavelength 0.9.*min_COEFF 1.1*max_COEFF]);
    
    subplot(1,5,[3,4,5]);
    imagesc(imrotate(SCORE(:,:,index),-angle));
    axis equal;
    axis off;
    colormap(gray)
    colorbar
    ax = gca;
    ax.FontSize = 20; 
    
    t2=title('Score map');
    t2.FontSize=16;
end

waitbar(2/2,hw2);
close (hw2);
    
    
    
end