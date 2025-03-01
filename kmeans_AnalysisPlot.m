function kmeans_AnalysisPlot(n_clusters,Mask_selection)
%This function plots the k-means clusters (number indicated by the user in
%n_clusters) and the correspondent centroids
%If Mask_selection is equal to 1 then the selection of a mask can be
%performed

MAT_HYP=0; %MAT_HYP=1 load .mat hypercube

if n_clusters>10
    n_clusters=10;
end

%colormap
% MAP_COLOR=[1 0 0;
%     0 1 0;
%     0 0 1;
%     1 0 1;
%     0 1 1;
%     1 1 0;
%     0.5 0 0.5;
%     0 0.5 0.5;
%     0.5 0.5 0;
%     0.3 0 0];

MAP_COLOR=[1 0 0;
    0 1 0;
    0 0 1;
    1 0 1;
    0 1 1;
    1 1 0;
    0.5 0 0.5;
    0 0.5 0.5;
    0.5 0.5 0;
    0.3 0 0];

if MAT_HYP==1
    dir0=('C:\Users\HARDi\Desktop\Tesi');
    [filename, pathname] = uigetfile('*.mat', 'Load Spectral Hypercube',dir0);
    file_tot=[pathname,filename];

    h=waitbar(1/2,'Loading the Spectra Hypercube');
    input_data=load(file_tot);
    close(h);
    Hyperspectrum_cube=abs(input_data.Hyperspectrum_cube);
    fr_real=input_data.fr_real;
    c=3*1e8;
    wavelength=c./fr_real.*1e-3;
    
    if isfield(input_data,'Intens')==1 %if Intens has been saved together with Hyperspectrum_cube **Benedetto** 12/07/2021
        Intens=input_data.Intens;
    else
        Intens=sum(abs(Hyperspectrum_cube),3);
    end
else
    [Hyperspectrum_cube,f,saturationMap,file_totCal,fr_real,Intens] = H5HypercubeRead;
    c=299792458; %speed of light [m/s]
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

%wavelength range selection
prompt={'Lower wavelength (nm)','Higher wavelength (nm)'};
name='Select spectral band and binning';
numlines=1;
defaultanswer={num2str(min_wavelength),num2str(max_wavelength),num2str(1)};

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer=inputdlg(prompt,name,numlines,defaultanswer,options);
lmCrop(1)=str2double(cell2mat(answer(1)));
lmCrop(2)=str2double(cell2mat(answer(2)));

f_Crop(1)=c./(lmCrop(1)*1e-9)/1e12; 
f_Crop(2)=c./(lmCrop(2)*1e-9)/1e12; 
f_Crop=sort(f_Crop);

Index_Crop=fr_real>f_Crop(1) & fr_real<f_Crop(2);

Hyperspectrum_cube=Hyperspectrum_cube(:,:,Index_Crop);
fr_real=fr_real(Index_Crop);
wavelength=c./fr_real.*1e-3;
max_wavelength=max(wavelength);
min_wavelength=min(wavelength);
%end wavelength range selection part

hw2=waitbar(1/2,'Hypercube k-means processing. Please wait...');

[CLUSTERS,CENTROIDS]=kmeans_tool(Hyperspectrum_cube,n_clusters);

max_CENTROIDS=max(max(CENTROIDS));
min_CENTROIDS=min(min(CENTROIDS));


pos1 = [100 50 1300 700];
pos2 = [0.1 0.3 0.3 0.3];

for index=1:n_clusters
    figure;
    set(gcf,'Position',pos1)
    
    subplot('Position',pos2)
    subplot(1,5,[1,2]);
    plot(wavelength,CENTROIDS(:,index),'Linewidth',2,'Color',MAP_COLOR(index,:));
    xlabel('wavelength [nm]', 'FontSize', 16);
    ylabel('amplitide', 'FontSize', 16);
    t1=title(['Centroid ',num2str(index)]);
    t1.FontSize=16;
    ax = gca;
    ax.FontSize = 16; 
    axis([min_wavelength max_wavelength 0.9.*min_CENTROIDS 1.1*max_CENTROIDS]);
    
    subplot(1,5,[3,4,5]);
    imagesc(imrotate(CLUSTERS==index,-angle));
    axis equal;
    axis off;
    MAP_index=[0 0 0; MAP_COLOR(index,:)]; 
    colormap(MAP_index)
    colorbar
    ax = gca;
    ax.FontSize = 20; 
    
    t2=title(['Cluster',num2str(index)]);
    t2.FontSize=16;
end

waitbar(2/2,hw2);
close (hw2);


%plot all the clusters together

figure;
set(gcf,'Position',pos1)

subplot('Position',pos2)
subplot(1,5,[1,2]);
for index=1:n_clusters
    plot(wavelength,CENTROIDS(:,index),'Linewidth',2,'Color',MAP_COLOR(index,:));
    hold on;
end
xlabel('wavelength [nm]', 'FontSize', 16);
ylabel('amplitide', 'FontSize', 16);
t1=title('Centroids');
t1.FontSize=16;
ax = gca;
ax.FontSize = 16; 
axis([min_wavelength max_wavelength 0.9.*min_CENTROIDS 1.1*max_CENTROIDS]);


subplot(1,5,[3,4,5]);
imagesc(imrotate(CLUSTERS,-angle));
axis equal;
axis off;
colormap(MAP_COLOR(1:n_clusters,:))
colorbar
ax = gca;
ax.FontSize = 20; 

t2=title('Clusters');
t2.FontSize=16;


%plot all clusters together with normalized centroids

figure;
set(gcf,'Position',pos1)

subplot('Position',pos2)
subplot(1,5,[1,2]);
for index=1:n_clusters
    plot(wavelength,CENTROIDS(:,index)./max(CENTROIDS(:,index)),'Linewidth',2,'Color',MAP_COLOR(index,:));
    hold on;
end
xlabel('wavelength [nm]', 'FontSize', 16);
ylabel('amplitide', 'FontSize', 16);
t1=title('Norm. Centroids');
t1.FontSize=16;
ax = gca;
ax.FontSize = 16; 
axis([min_wavelength max_wavelength 0 1]);


subplot(1,5,[3,4,5]);
imagesc(imrotate(CLUSTERS,-angle));
axis equal;
axis off;
colormap(MAP_COLOR(1:n_clusters,:))
colorbar
ax = gca;
ax.FontSize = 20; 

t2=title('Clusters');
t2.FontSize=16;
    
    
    
end