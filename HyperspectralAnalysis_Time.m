% Cristian Manzoni 4.2.2019, from Beam Profiles
% Benedetto Ardini 23.06.2023
% Martina Riva 07.04.2025



clear all


% --VERSIONE LINUX O WINDOWS-- 1
if isunix
    linux=1; %0 se si usa Windows, 1 (o altri valori) se si usa Linux
else
    linux=0;
end;

if ne(linux,0) %se si sta usando Linux
    fprintf('\n\n--LINUX VERSION of HyperspectralAnalysis_Spectrumfi--\n\n');
else %se non si sta usando linux
    fprintf('\n\n--WINDOWS VERSION of HyperspectralAnalysis_Spectrum--\n\n');
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
    path(path,'C:\Users\marti\OneDrive - Politecnico di Milano\PhD\Measurements');

    dir0=('C:\Users\marti\OneDrive - Politecnico di Milano\PhD\Measurements');
end

%**Benedetto**

home;

fprintf(['\n\n *** Analysis of Hyperspectral image ***']);

load_matrix=round(Dinput('\n\n    Load pre-processed hypercube [0: no - 1: yes]? ',0));

if load_matrix

    [filename_Hyper, pathname_Hyper] = uigetfile('*.mat;*.h5;', 'Load Hypercube',dir0);
    file_tot=[pathname_Hyper,filename_Hyper];
    dir2=pathname_Hyper;

    stringa=[' *** Analysis of file ',file_tot,' ***'];
    fprintf('\n\n%s',stringa);

    if strcmp(filename_Hyper(end-2:end),'.h5') %for .h5 file
      H5toMatlabHypercube(filename_Hyper, pathname_Hyper); %creates a .mat in the same folder 
      file_tot=[pathname_Hyper, strcat(filename_Hyper(1:end-3), '_hyp.mat')]; %file to be loaded
    end

    load(file_tot);

    HyperMatrix=double(HyperMatrix); %HyperMatrix could have been saved in single

    % t=t(:,1)';
    pnt=size(HyperMatrix,3);
    fprintf(['\n\n  -- Number of images: ',num2str(pnt)]);

    step=max(diff(t));


    Averages=mean(HyperMatrix,3);
    Sample=abs(Averages)/(max(max(abs(Averages)))); %era senza abs e c'era .*256 **Benedetto** 30/06/2021

    A=HyperMatrix(:,:,1);
    [a,b]=size(A);
    fprintf(['\n  -- Size: ',num2str(a),' x ',num2str(b),' pixels']);


    %%%%%%**Benedetto** 14/09/2020
    % Linear scale
    hh=figure;
    %Sample_image is better to be visualized with the image plot **Benedetto** 30/06/2021
    Sample_image(:,:,1)=Sample;
    Sample_image(:,:,2)=Sample;
    Sample_image(:,:,3)=Sample;
    image(Sample_image); colormap(MAP); %era Sample (l'immagine di Sample veniva troppo scura) **Benedetto** 30/06/2021
    clear Sample_image

    title('Zoom and press ENTER. Only zoomed area will be transformed');
    pause;
    V=round(axis);

    close(hh);

    rmin=round(V(3)); rmax=round(V(4))-1; % First and last row
    cmin=round(V(1)); cmax=round(V(2))-1; % First and last column

    A=A(rmin:rmax,cmin:cmax);
    [a,b]=size(A);
    HyperMatrix=HyperMatrix(rmin:rmax,cmin:cmax,:); %the HyperMatrix image will be reduced accordin to the Zoom
    Averages=mean(HyperMatrix,3);
    Sample=abs(Averages)/(max(max(abs(Averages))))*256; %era senza abs

    fprintf(['\n  -- Number of selected pixels: ',num2str(a),' x ',num2str(b),' = ',num2str(a*b),' --']);

    %%%%%%**Benedetto** 14/09/2020


    %%%%%%%%%%%%%%%%%%
    Corr_summ=zeros(pnt,100);
    counter=0;

    y_samp=round(linspace(1,a,10));
    x_samp=round(linspace(1,b,10));

    for y=1:10

        for x=1:10
            counter=counter+1;
            Corr_summ(:,counter)=squeeze(HyperMatrix(y_samp(y),x_samp(x),:)-Averages(y_samp(y),x_samp(x)));
            % Corr_summ(:,counter)=Corr_summ(:,counter)/max(Corr_summ(:,counter)); % Normalizzata

        end;

    end;

    hh=figure;
    plot(Corr_summ,'linewidth',2); axis tight;
    title('Zoom and press ENTER. Only zoomed area will be transformed');
    zoom xon;
    pause;
    V=round(axis);



    pnt=V(2)-V(1)+1; % New value of pnt

    fprintf(['\n  -- A subset of ',num2str(V(2)-V(1)+1),' images will be saved --']);

    %%%%%%**Benedetto** 15/09/2020

    pola=Dinput('\n\n    Invert interferogram amplitude: (0 = no, 1 = yes): ',1);
    close(hh);

    if pola % if polarizers are crossed, data should be flipped
        HyperMatrixSwap=-HyperMatrix;
        HyperMatrix=HyperMatrixSwap;
        clear HyperMatrixSwap;
    end;

    %%%%%%**Benedetto** 15/09/2020

    % Removing unnecessary data
    HyperMatrix=HyperMatrix(:,:,V(1):V(2));
    t=t(V(1):V(2));

    Averages=mean(HyperMatrix,3);
    Sample=abs(Averages)/(max(max(abs(Averages))))*256; %era senza abs

else

    dir2=uigetdir(dir0,'Select Directory');

    form=Dinput('\n File format (0: TIFF - 1: binary Milan - 2: binary BH)',1);

    el=dir(dir2);
    fine=size(el,1)-2; % number of files in folder

    ii=1;
    while (el(ii).isdir)
        ii=ii+1;
    end;

    stringa=[' *** Analysis of file ',el(ii).name,' ***'];
    fprintf('\n\n%s',stringa);

    switch  form
        case 0
            A=double(([dir2,'/',el(ii).name]))/256;
        case 1
            fprintf('\n\n Sensor size:');
            fprintf(  '\n   1. 1004 x 1002 \n   2. 1024 x 1280 \n   3. 1280 x 1024 \n   4. Other');
            Sensor=Dinput('\n\n Select sensor size:',1);

            switch Sensor
                case 1
                    Sx=1004;
                    Sy=1002;
                case 2
                    Sx=1024;
                    Sy=1280;
                case 3
                    Sy=1024;
                    Sx=1280;
                case 4
                    Sx=Dinput('\n X dimension:',1004);
                    Sy=Dinput('\n Y dimension:',1002);
            end

            fileID = fopen([dir2,'/',el(ii).name]);
            A = fread(fileID,[Sx Sy],'uint32');
            fclose(fileID);
        case 2
            A=load([dir2,'/',el(ii).name]);
    end;

    A=A(:,:,1); % Takes the first useful file.

    hcut=figure;
    %     % log scale
    %     yy=log((A-min(min(A))));
    %     yy=yy/max(max(yy))*256;
    %     image(yy); colormap(MAP); % axis equal;

    % linear scale
    imagesc(A); colormap(MAP); % axis equal;
    axis tight; axis off;

    % % Linear scale
    % image(A); colormap(MAP); axis equal; axis off;

    title('Zoom and press ENTER. Only zoomed area will be transformed');
    pause;
    V=round(axis);

    rmin=round(V(3)); rmax=round(V(4))-1; % First and last row
    cmin=round(V(1)); cmax=round(V(2))-1; % First and last column

    A=A(rmin:rmax,cmin:cmax);
    [a,b]=size(A);

    fprintf(['\n  -- Number of selected pixels: ',num2str(a),' x ',num2str(b),' = ',num2str(a*b),' --']);

    Bin=Dinput('\n\n    Binning NxN. Select N (1: no binning)',1);
    Bin=ceil(abs(Bin));
    aBin=floor(a/Bin);
    bBin=floor(b/Bin);

    fprintf(['\n  -- Number of produced pixels: ',num2str(aBin),' x ',num2str(bBin),' = ',num2str(aBin*bBin),' --']);

    %     % Generation of the noise matrix: loading file from folder /Bkg
    %     noise_removal=round(Dinput('\n\n    Remove noise? [0: No - 1: From File - 2: constant] ',0));
    %
    %     switch noise_removal
    %
    %         case 1
    %
    %             dirBkg=[dir2,'/Bkg'];
    %             elBkg=dir(dirBkg);
    %
    %             ii=1;
    %             while (elBkg(ii).isdir)
    %                 ii=ii+1;
    %             end;
    %
    %             switch  form
    %                 case 0
    %                     B=double(imread([dirBkg,'/',elBkg(ii).name]))/256;
    %                 case 1
    %                     fileID = fopen([dirBkg,'/',elBkg(ii).name]);
    %                     B = fread(fileID,[Sx Sy],'uint32');
    %                     fclose(fileID);
    %                 case 2
    %                     B=load([dirBkg,'/',elBkg(ii).name]);
    %             end;
    %
    %             noise=B(rmin:rmax,cmin:cmax,1); % Takes the background.
    %
    %         case 2
    %
    %             n_bkg=(Dinput('\n   Noise value: ',0));
    %             noise=ones(size(A))*n_bkg;
    %
    %         case 0
    %
    %             noise=zeros(size(A));
    %
    %     end;

    %%%%%%%%%%%
    % BINNING: theResult = sepblockfun(A(1:aBin*Bin,1:bBin*Bin),[Bin,Bin],'mean');
    %%%%%%%%%%%

    %    noiseN = sepblockfun(noise(1:aBin*Bin,1:bBin*Bin),[Bin,Bin],'mean');

    pnt=0;

    %    HyperMatrix=zeros([size(noiseN),fine]); % Smaller matrix

    h=waitbar(0,'Loading all files');

    switch  form
        case 0
            for ii=1:size(el,1)

                if not(el(ii).isdir)

                    pnt=pnt+1;

                    A=double(imread([dir2,'/',el(ii).name]))/256;
                    A=A(rmin:rmax,cmin:cmax,1); % Takes the first matrix.

                    % A=A-noise; % Removes the noise
                    % HyperMatrix(:,:,pnt)=A;
                    HyperMatrix(:,:,pnt) = sepblockfun(A(1:aBin*Bin,1:bBin*Bin),[Bin,Bin],'mean');

                    waitbar(pnt/fine,h);

                end;

            end
        case 1

            for ii=1:size(el,1)

                if not(el(ii).isdir)

                    pnt=pnt+1;

                    fileID = fopen([dir2,'/',el(ii).name]);
                    A = fread(fileID,[Sx Sy],'uint32');
                    fclose(fileID);

                    A=A(rmin:rmax,cmin:cmax,1); % Takes the first matrix.

                    % A=A-noise; % Removes the noise
                    % HyperMatrix(:,:,pnt)=A;
                    HyperMatrix(:,:,pnt) = sepblockfun(A(1:aBin*Bin,1:bBin*Bin),[Bin,Bin],'mean');

                    waitbar(pnt/fine,h);

                end;

            end

        case 2

            for ii=1:size(el,1)

                if not(el(ii).isdir)

                    pnt=pnt+1;

                    A=load([dir2,'/',el(ii).name]);

                    A=A(rmin:rmax,cmin:cmax,1); % Takes the first matrix.

                    % A=A-noise; % Removes the noise
                    % HyperMatrix(:,:,pnt)=A;
                    HyperMatrix(:,:,pnt) = sepblockfun(A(1:aBin*Bin,1:bBin*Bin),[Bin,Bin],'mean');

                    waitbar(pnt/fine,h);

                end;

            end
    end;

    fprintf(['\n  -- Number of images: ',num2str(pnt),' --']);

    a=aBin; b=bBin; % New matrix dimension, after binning

    waitbar(1,h,'Preprocessing');

    % pnt: size of each correlation
    HyperMatrix=HyperMatrix(:,:,1:pnt);

    Averages=mean(HyperMatrix,3);
    Averages0=Averages;
    Sample=abs(Averages)/(max(max(abs(Averages))))*256; % Will be used to show the image %era senza abs
    close(h);

    % Generation of the time axis: loading file from folder
    time_axis=round(Dinput('\n\n    Load time axis [0: no - 1: yes]? ',0));

    if time_axis

        dirTime=uigetdir(dir2,'Select Directory');
        elTime=dir(dirTime);

        ii=1;
        while (elTime(ii).isdir)
            ii=ii+1;
        end;

        t=load([dirTime,'/',elTime(ii).name]);

        t=t(1:pnt); % <-- sistemare l'acquisizione
        if iscolumn(t), t=t'; end;

        step=max(diff(t));
        fprintf(['\n  -- Average sampling step: ',num2str(mean(diff(t))),' --']);

        units=round(Dinput('\n    Units [3: mm - 6: micron]? ',3));
        t=t.*10^(6-units); % Converted into microns

        step=max(diff(t));
        fprintf(['\n  -- Average sampling step: ',num2str(mean(diff(t))),' micron --']);

    else

        step=Dinput('\n    Sampling step: ',1');

        % Generating time axis
        % t=[1:pnt]-(pnt-1)/2;
        t=([1:pnt]-(pnt-1)/2)*step;
    end;

    % Corr_summ=zeros(pnt,round(a*b/100));
    Corr_summ=zeros(pnt,100);

    % Per eliminare alcuni elementi non necessari

    counter=0;

    y_samp=round(linspace(1,a,10));
    x_samp=round(linspace(1,b,10));

    for y=1:10

        for x=1:10
            counter=counter+1;
            Corr_summ(:,counter)=squeeze(HyperMatrix(y_samp(y),x_samp(x),:)-Averages(y_samp(y),x_samp(x)));
            % Corr_summ(:,counter)=Corr_summ(:,counter)/max(Corr_summ(:,counter)); % Normalizzata

        end;

    end;

    hh=figure;
    plot(Corr_summ,'linewidth',2); axis tight;
    title('Zoom and press ENTER. Only zoomed area will be transformed');
    zoom xon
    pause;
    V=round(axis);


    %     Corr_summ=Corr_summ(V(1):V(2),:);
    %     plot(Corr_summ,'linewidth',2); axis tight;
    %     title('Selected correlations. Press ENTER');
    %     pause;
    close(hh);

    pnt=V(2)-V(1)+1; % New value of pnt

    fprintf(['\n  -- A subset of ',num2str(V(2)-V(1)+1),' images will be saved --']);

    % Removing unnecessary data
    HyperMatrix=HyperMatrix(:,:,V(1):V(2));
    t=t(V(1):V(2));
    fprintf(['interferogram extremes are',num2str(V(1)),' ',num2str(V(2))]);
    %%%%%%%%%%%%%% Era sotto

    pola=Dinput('\n\n    Invert interferogram amplitude: (0 = no, 1 = yes): ',1);

    if pola % if polarizers are crossed, data should be flipped
        HyperMatrixSwap=-HyperMatrix;
        HyperMatrix=HyperMatrixSwap;
        clear HyperMatrixSwap;
    end;

    reverse=Dinput('    Reverse Delay Axis: (0 = no, 1 = yes): ',0);

    if reverse

        h=waitbar(0,'Reversing Time axis');
        HyperMatrixSwap=zeros(size(HyperMatrix));

        for nn=1:pnt

            waitbar(nn/(pnt+10),h);
            HyperMatrixSwap(:,:,nn)=HyperMatrix(:,:,pnt-nn+1);
        end;

        HyperMatrix=HyperMatrixSwap;
        clear HyperMatrixSwap;

        close(h);

    end;

    %%%%%%%%%%%%%%%
    %Start of the code SPATIAL SHIFT TEMPORAL HYPERCUBE CORRECTION **Benedetto** 27/05/2022
    %Code added for spatial shift correction along y-axis of the image
    %Updated with gradual correction on 29/06/2022
    fprintf('\n\n    Spatial y-axis shift correction (necessary for YVO4 interferometer):');
    fprintf(  '\n    1. standard YVO4+LucaR camera correction \n    2. custom correction \n    3. no correction');
    SpatialShiftCorrection=Dinput('\n\n    Select choice:',1);

    if SpatialShiftCorrection==1 || SpatialShiftCorrection==2

        if SpatialShiftCorrection==1
            wedge_excursion=17999; %YVO4 wedge excursion [um] in the 220527\USAF_0 measurement with LucaR camera
            spatial_shift=17; %spatial shift [pixels] for the wedge excursion above
            pixel_ratio=1002./Sy.*Bin; %image pixel ratio with respect to the LucaR 8um pixel
            %(1002/Sy): binning in detection (pixel dimension compared to the LucaR pixel)
            % Bin: software binning
        else
            wedge_excursion=Dinput('\n\n    Wedge excursion [um]:',9600);
            spatial_shift=Dinput('\n\n    Related image spatial shift [pixels] (without detection binning):',5);
            %pixel_dimension=Dinput('\n\n    Camera pixel dimension [um]:',8);
            binning=Dinput('\n\n    Camera binning applied in detection:',1);
            pixel_ratio=binning.*Bin;
            %pixel_ratio=pixel_dimension./8.*binning.*Bin; %image pixel ratio with respect to the LucaR 8um pixel
            %clear pixel_dimension binning
        end

        shift_rate=spatial_shift./wedge_excursion./pixel_ratio; %pixels shift velocity [pixels/um]
        shift_num=floor((max(t)-min(t)).*shift_rate); %number of cumulative 1-pixel shift to apply during the delay scan

        transl=0; %pixels translation
        R_pxl_w=0; %right (adjacent) pixel weight for gradual shift correction
        L_pxl_w=1; %left (current) pixel weigth for gradual shift correction

        h_shiftCorrection=waitbar(0,'Temporal hypercube spatial shift correction...');

        %gradual shift correction of temporal hypercube frames
        for index=2:length(t) %the first frame is not shifted (index start from 2)

            if R_pxl_w>=1 %once the left (adjacent) pixel weight is greater than one it means that you need to apply a translation 1 pixel on the image
                R_pxl_w=R_pxl_w-1;
                L_pxl_w=1-R_pxl_w;
                transl=transl+1; %update transl by adding 1 (you need a translation of 1 pixel more)
            end

            R_pxl_w=R_pxl_w+shift_rate.*(t(index)-t(index-1)); %the right (adjacent) pixel weight is calculated considering sub-pixel image shift due to the wedge step t(index)-t(index-1) with respect to the pixel size
            L_pxl_w=1-R_pxl_w;

            %correction of the frame n° index by a weighted average of the
            %frame (matrix of current pixels) with its copy shifted by 1
            %pixel on the y direction (matrix of adjacent pixels)
            HyperMatrix(:,1:end-shift_num-1,index)=R_pxl_w.*HyperMatrix(:,2+transl:end-shift_num+transl,index)+L_pxl_w.*HyperMatrix(:,1+transl:end-shift_num-1+transl,index);

            waitbar((index-1)./(length(t)-1),h_shiftCorrection);

        end

        HyperMatrix=HyperMatrix(:,1:end-shift_num-1,:); %pixels in the border are discarded
        b=size(HyperMatrix,2); %new y dimension of the matrix
        close(h_shiftCorrection);
        spatialShiftCorrection=['Spatial shift correction on Y axis applied, shift rate = shift/wedge excursion =' num2str(shift_rate) ' [pixels/um]']; %string to save with Temporal Hypercube
        clear h_shiftCorrection time_interval shift_rate shift_num

    else
        spatialShiftCorrection='No spatial shift correction applied'; %string to save with Temporal Hypercube
    end
    %End of the code SPATIAL SHIFT TEMPORAL HYPERCUBE CORRECTION
    %%%%%%%%%%%%%%%

    Averages=mean(HyperMatrix,3);
    Averages0=Averages;  % Aggiunto

    % Saving image and step
    [filename_Hyper, pathname_Hyper] = uiputfile('*.mat', 'Save Hypercube as',dir2);
    file_tot=[pathname_Hyper,filename_Hyper];

    stringa=['  -- Saving Hypercube in ',file_tot,' --'];
    fprintf('\n%s',stringa);

    % Checks the size of the matrix, before saving
    HyProp = whos('HyperMatrix') ;
    Gigabytes = HyProp.bytes/2^30;
    settings=struct('Amplitude_inversion', pola, 'Time_axis_inversion', reverse, 'Axis_limits', [V(1), V(2)], 'Temporal_axis', t);

    if Gigabytes<1.8
        save(file_tot,'t','HyperMatrix','Averages0','spatialShiftCorrection'); % raw data, averages not removed yet
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Code added for temporal hypercube settings saving- Martina 241223
        save(strcat(file_tot, '_settings'), 'settings'); %settings in temp hyp processing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        save(file_tot,'t','HyperMatrix','Averages0','spatialShiftCorrection','-v7.3'); % raw data, averages not removed yet
        save(strcat(file_tot, '_settings'), 'settings'); %settings in temp hyp processing
    end;

    % save(file_tot,'t','HyperMatrix','Averages0','-v7.3'); % raw data, averages not removed yet

end;

% End of preparation of hypercube


%%%**Benedetto** CREATION OF MAP OF SATURATED PIXELS (0:saturated, 1:not saturated)
%%%modifica: 25/10/2020
fprintf('\n\n Removing saturated pixels. Please select sensor dynamics:');
fprintf(  '\n   1. 14 bit \n   2. 12 bit \n   3. Other');
dynamics=Dinput('\n\n Select sensor dynamics:',1);
switch dynamics
    case 1
        saturationThreshold=2.^(14).*(1-0.05); %put the threshold 5% levels below the true saturation level
    case 2
        saturationThreshold=2.^(12).*(1-0.05); %put the threshold 5% levels below the true saturation level
    case 3
        n_bit_dynamics=Dinput('\n\n Number of bits of sensor dynamics:',14);
        saturationThreshold=2.^(n_bit_dynamics).*(1-0.05); %put the threshold 5% levels below the true saturation level
        clear n_bit_dynamics;
end

saturationMap=ones(size(HyperMatrix,1),size(HyperMatrix,2));
maxHyperMatrix_dynamics=max(abs(HyperMatrix),[],3); %maximum calculated on the time dimension
saturationMap(maxHyperMatrix_dynamics>saturationThreshold)=0; %saturated pixels are mapped with 0s

Sample_dynamics=abs(Averages)./max(max(abs(Averages))); %era senza abs

hh1=figure;
Pixels=zeros(size(Sample_dynamics,1),size(Sample_dynamics,2),3);
Pixels_1=Sample_dynamics.*255;
Pixels_2=Sample_dynamics.*255;
Pixels_3=Sample_dynamics.*255;
Pixels_1(saturationMap==0)=255;
Pixels_2(saturationMap==0)=0;
Pixels_3(saturationMap==0)=0;
Pixels(:,:,1)=Pixels_1;
Pixels(:,:,2)=Pixels_2;
Pixels(:,:,3)=Pixels_3;
hhh1=image(uint8(Pixels));
title('Saturated pixels are shown in red (click ENTER to close)');
pause;
close(hh1);
clear maxHyperMatrix_dynamics Sample_dynamics;
%clear hh1 hhh1 Pixels_1 Pixels_2 Pixels_3 Pixels; %original version

% **Martina**-11.02.2025: save map of saturated pixels
disp('Saving Saturated pixels map')
save(strcat(pathname_Hyper, 'SaturatedPixelsMap'), 'saturationMap');
clear hh1 hhh1 Pixels_1 Pixels_2 Pixels_3 Pixels;
%%**Benedetto** END of MAP OF SATURATED PIXELS section


% Removal of averages
for nn=1:pnt
    HyperMatrix(:,:,nn)=HyperMatrix(:,:,nn)-Averages;
end;


% Extraction of a subset of interferograms

Corr_summ=zeros(pnt,100);
counter=0;

y_samp=round(linspace(1,a,10));
x_samp=round(linspace(1,b,10));

for y=1:10

    for x=1:10
        counter=counter+1;
        if y_samp(y)<size(HyperMatrix,1) && y_samp(x)<size(HyperMatrix,2) %if-else **Benedetto 27/05/2022
            Corr_summ(:,counter)=squeeze(HyperMatrix(y_samp(y),x_samp(x),:)); %-Averages(y_samp(y),x_samp(x)));
        else
            Corr_summ(:,counter)=squeeze(HyperMatrix(size(HyperMatrix,1),size(HyperMatrix,2),:));
        end

    end;

end;

%
% Gif generation

if Dinput('\n\n    Generate time GIF (0 = No, 1 = Yes)? ',0)

    n_frames=Dinput('\n    Number of frames: ',100);

    Maximum=max(max(max(HyperMatrix)));

    [filename_GIF, pathname_GIF] = uiputfile('*.gif', 'Save GIF as',dir0);
    filename=[pathname_GIF,filename_GIF];

    n_max=length(t);
    nn=round(linspace(1,n_max,100));

    figure(1)
    for n = 1:n_frames

        image(HyperMatrix(:,:,nn(n))/Maximum*256*10); axis equal; axis off;
        rectangle('FaceColor',[1 0 0 ],'EdgeColor','none','Position',[b*min([n n_frames/2])/n_frames 0 abs(b/n_frames*(n-n_frames/2)) a/30 ]);
        colormap(MAP);
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if nn(n) == 1;
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end

end;
%%%%%%%%%%%%%%

% Amplitude=max(abs(HyperMatrix(:,:,1:pnt)),[],3);
% Normalization=Averages./Amplitude; % Required to normalize the spectrum to the amplitude

h0=figure(1);
subplot('position',[0 0.05 0.7 0.9]);
h0a=gca; axis off;

subplot('position',[0.72 0.5 0.28 0.45]);
h0b=gca; axis off;

subplot('position',[0.72 0.05 0.28 0.4]);
h0c=gca; axis off;

fprintf('\n\n    Delay correction:'); %**Benedetto** 29/01/2021
fprintf(  '\n    1. load delay correction file \n    2. select spectra for delay correction \n    3. no delay correction');
DelayCorrection=Dinput('\n\n    Select choice:',1);

if DelayCorrection==1 %load external calibration for delay correction **Benedetto** 29/01/2021
    [filename_delayCal, pathname_delayCal] = uigetfile('*.mat', 'Load Delay Correction',dir0);
    file_totDelayCal=[pathname_delayCal,filename_delayCal];

    stringa=[' *** Loaded Delay Correction ',file_totDelayCal,' ***'];
    fprintf('\n\n%s',stringa);

    load(file_totDelayCal); % Contains t_axis and  t_delayCorrection

    t_correzione=interp1(t_axis,t_delayCorrection,t); %the delay correction is referred to t_axis: it has to be interpolated in the t axis of the specific measurement
    t_correzione(isnan(t_correzione))=0;

    tt=t+t_correzione;

    h0=figure(1); axes(h0a);

    % % linear scale
    imagesc(Sample.^0.5); axis equal; axis off
    colormap(gray(256));
else
    if DelayCorrection==2 %delay correction with selected spectra %**Benedetto** 24/06/2021
        sp_corr=round(Dinput('\n    How many spectra for delay correction (0 for no correction)? ',0));
    else %no delay correction (delay correction with 0 spectra) %**Benedetto** 24/06/2021
        sp_corr=0;
    end

    if sp_corr>0

        Ns=floor(Dinput('\n    Delay correction: Binning over N x N neighboring pixels: N = ', 1)/2);
        % We will average in the range x-Ns:x+Ns , y-Ns:y+Ns  <- Check

        t_corretto=zeros(sp_corr,pnt);
        t_correzione=zeros(sp_corr,pnt);

        conversionWindow=WConversion_cryst; %finestra di aiuto per le conversioni **BENEDETTO**
        set(conversionWindow,'OuterPosition', [205 -7 114.5714 35.1765]);

        fprintf('\n\n -- Select the frequency range and resolution for spectrum--');
        fmin=Dinput('\n\n    Minimum pseudo-frequency: ',0);
        fmax=Dinput('    Maximum pseudo-frequency: ',1/(2*step));
        resol=round(Dinput('    Number of spectral points: ',500));

        close(conversionWindow);
        %f=linspace(0,1/(2*step),501); %prima era questo **Benedetto**
        %28/01/2021
        f=linspace(fmin,fmax,resol);

        figure(h0); axes(h0a);

        %  linear scale
        % image(Amplitude.^0.5); axis equal; axis off;
        imagesc(Sample.^0.5); axis equal; axis off;

        colormap(MAP);

        % ha=gca;
        title('Linear amplitude - Gamma 0.5 - Zoom if necessary, then ENTER');
        pause;

        for hh=1:sp_corr

            figure(h0); axes(h0a);
            % pause;
            [x,y]=ginput(1);

            x=round(x);
            y=round(y);

            % C=squeeze(HyperMatrix(y,x,:))';
            C=squeeze(mean(mean((HyperMatrix(y-Ns:y+Ns,x-Ns:x+Ns,:)),2),1))';

            C1=C-mean(C);
            filtro=Apodization(4,length(C1),length(C1)/2); % Era 1
            C10=C1.*filtro;

            figure(h0); axes(h0b);
            plot(t,C1,'r',t,C10,'b','linewidth',2);

            SpettroC10=FourierDir(t,C10,f);

            figure(h0); axes(h0c);
            plot(f,abs(SpettroC10),'b','linewidth',2);
            legend('Zoom and press ENTER. Only zoomed area will be back-transformed');
            pause;
            V=axis;

            SpettroC20=SpettroC10;

            f11=linspace(V(1),V(2),1001); % Asse delle frequenze limitato alla regione zoomata
            SpettroC10=FourierDir(t,C10,f11); % *  Trasformo solo la parte che mi serve. Eventualmente questa riga si puÃ² togliere, togliere quelle con * e rimettere quelle con **
            % SpettroC10(f<V(1))=0; SpettroC10(f>V(2))=0;  **

            [~,xp_coord]=max(abs(SpettroC10));
            % xp=f(xp_coord); **
            xp=f11(xp_coord); % *

            % Ricavo la correzione dell'asse dei tempi
            % yt=FourierInv(f,SpettroC10,t); % **
            yt=FourierInv(f11,SpettroC10,t); % *

            phi=unwrap(angle(yt));
            tt=phi./(2*pi*xp); tt=tt-tt(1)+t(1); % Nuovo asse
            t_corretto(hh,:)=tt;
            t_correzione(hh,:)=tt-t;

            C30=C1.*filtro;

            % New spwctrum after correction
            % SpettroC30=FourierDir(tt,C30,f); **
            SpettroC30=FourierDir(tt,C30,f11); % *

            figure(h0); axes(h0c);
            % plot(f,abs(SpettroC20),'b',f,abs(SpettroC30),'r','linewidth',2); **
            plot(f11,abs(SpettroC10),'b',f11,abs(SpettroC30),'r','linewidth',2); % *
            legend('Original','Corrected');

        end;

        tt=mean(t_corretto,1);
        tt_plot=mean(t_correzione,1);

        figure(h0); axes(h0a);
        plot(t,t_correzione',t,tt_plot,'k--','linewidth',2);
        legend('Time axis correction');

        delayCorrectionSave=menu('Save delay correction?','Yes','No'); %**Benedetto 28/01/2021
        if delayCorrectionSave==1 %save t, the delay axis, and tt_plot, the mean of the delay corrections in t_correzione
            [filename_DelayCorr, pathname_DelayCorr] = uiputfile('*.mat', 'Save Delay Correction',dir0);
            file_DelayCorr=[pathname_DelayCorr,filename_DelayCorr];
            save(file_DelayCorr,'t','tt_plot');
        end

    else
        fprintf('\n    No delay correction will be applied\n');

        tt=t;
        h0=figure(1); axes(h0a);

        % % linear scale
        imagesc(Sample.^0.5); axis equal; axis off
        colormap(gray(256));

    end
end

preview=Dinput('\n\n Do you want to see a preview for selected pixels spectra (Yes [1]; No [0])? ', 0); 
% there are 2 possibilities:
%     [1]-> you can load the calibration or make it 
%     [0]-> you avoid preview, but you must have calibration

switch preview
    case 1 %Preview

        fprintf('\n  -- Spectra of preview pixels --');

        conversionWindow=WConversion_cryst; %finestra di aiuto per le conversioni **BENEDETTO**
        set(conversionWindow,'OuterPosition', [205 -7 114.5714 35.1765]);

        % sp=round(Dinput('\n\n    How many preview pixels? ',3));
        fmin=Dinput('\n\n    Minimum pseudo-frequency: ',0);
        fmax=Dinput('    Maximum pseudo-frequency: ',1/(2*step));
        resol=round(Dinput('    Number of spectral points: ',500));

        close(conversionWindow);

        linearFit=Dinput('\n    Subtraction from interferogram of: [0: constant line - 1: linear trend- 2: quadratic trend]? ',0);
        trapz=round(Dinput('\n    Apodization (0: none - 1: Happ-Genzel - 2: 3-term Blackman-Harris - 3: 4-term Blackman-Harris - 4: Triangular - 5: Supergaussian - 6: Trapezoidal)? ',1)); %**Benedetto**
        center_find=round(Dinput('\n    Interferograms center calculation (0: exactly in the middle   - 1: function baricenter)? ',1)); %**Benedetto** 14/01/2021

        prompt={'Pixels binning for the mean'}; %Possibility to select interferogram mean over more than one pixel **Benedetto** 11/01/2022
        name='Select pixels mean';
        numlines=1;
        defaultanswer={num2str(1)};

        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';

        answer=inputdlg(prompt,name,numlines,defaultanswer,options);

        bin_preview=str2double(cell2mat(answer));
        if bin_preview==0
            bin_preview=1;
        end

        %%%**Benedetto** ha messo come commento le seguenti righe il 09/10/2020
        % % %**Benedetto** 26/09/2020
        % % realPart=Dinput('\n    SPECTRUM: [0: absolute value (asymmetric interf., e.g. Raman) - 1: real part (symmetric interf.)]? ',1);
        % % %**Benedetto** 26/09/2020

        if trapz==6 %**Benedetto**

            auto_max=Dinput('\n    Trapezoidal apodization. Retrieval of maximum: [0: automatic - 1: manual]',0);

            if auto_max==1
                hh=figure;

                Aver_correl=sum(Corr_summ,2); % Average correlation

                % plot(Corr_summ,'linewidth',2); axis tight;
                plot(Aver_correl,'linewidth',2); axis tight;

                zoom xon; pause;
                [Center,~]=ginput(1);

                close(hh);

            end;

        end;

        % ref=round(Dinput('\n    The first point will be a reference [0: no - 1: yes]? ',0));

        tmax=tt(end);
        tmaxindex=find(tt<=tmax,1,'last'); % Index of the last element

        figure(h0); axes(h0a);
        % % linear scale
        imagesc(Sample.^0.5); axis equal; axis off % it was Sample
        colormap(gray(256));

        title(sprintf('Gamma: 0.5 - Zoom if necessary, then ENTER - RIGHT click when finished\nMean over %d x %d pixels selected',bin_preview,bin_preview));

        f=linspace(fmin,fmax,resol);

        tt0=tt;
        tt=tt0(1:tmaxindex);

        Dt=diff(tt);
        Dt(end+1)=Dt(end);
        exps=exp(-1i*2*pi*tt'*f);

        % for hh=1:sp

        hh=0;

        figure(h0); axes(h0a);
        pause; [x,y,button]=ginput(1);

        while button==1

            hh=hh+1;
            x=round(x);
            y=round(y);


            if mod(bin_preview,2)==0 %bin_preview even
                C0=squeeze(mean(HyperMatrix(y-(bin_preview/2-1):y+bin_preview/2,x-(bin_preview/2-1):x+bin_preview/2,:),[1,2]))';
            else %bin_preview odd
                C0=squeeze(mean(HyperMatrix(y-floor(bin_preview/2):y+floor(bin_preview/2),x-floor(bin_preview/2):x+floor(bin_preview/2),:),[1,2]))';
            end

            % cut of the trace
            % C0=C0(1:tmaxindex);
            C=C0;
            % C=medfilt1(C0,3); % removes glitch
            % C=medfilt1(C,3); % once again

            if linearFit  % Subtracts a line instead of a constant

                FitPar = polyfit(tt,C,1);

                FitLinear = polyval(FitPar,tt);
                C=C-FitLinear;

            elseif linearFit==2

                FitPar = polyfit(tt,C,2);

                FitQuadr = polyval(FitPar,tt);
                C=C-FitQuadr;
            end

            %**Benedetto**
            if trapz==6 %the CASE OF TRAPEZOIDAL apodization is treated independently **Benedetto**

                if auto_max==1

                    I=round(Center);
                    [maxC,~] = max(C);
                else

                    [maxC,I] = max(C); % finds the position of the max

                end;

                I1=find(tt-tt(I)<=(tt(I)-tt(1))); % tt(I)-tt(1) is the distance between the first
                % point and the peak. It will be the distance
                % between the peak and the last point

                tt1=tt(I1); % row
                y1=C(I1); % I select only the symmetric correlation

                minC = min(C); % finds the min

                maxC=maxC+Averages(y,x);
                minC=minC+Averages(y,x);

                Contrast(hh)=(maxC-minC)/(maxC+minC);
                Ave(hh)=Averages(y,x);
                Correl(hh,:)=(1-2*pola)*C+Averages(y,x); % Keeps the correlation; (1-2*pola) is -1 or +1, depending whether the correlation has been reversed

                Dt_1=diff(tt1);
                Dt_1(end+1)=Dt_1(end);

                exps_1=exp(-1i*2*pi*tt1'*f);

                s_1=y1.*Apodization(4,length(y1),I); % row

                % fourier computed with the exponential matrix and no loop involved
                yf_1=(Dt_1.*s_1)*exps_1;

                %%%**Benedetto** ha messo come commento le seguenti righe il 09/10/2020
                % %             phase=angle(yf_1); % This phase will be used to correct the spectral phase of the entire trace
                % %
                % %             if realPart %**Benedetto** 26/09/2020
                % %                 %Interpolazione mediante least mean squares
                % %                 xxx=f;
                % %                 Amp=abs(yf_1).^2;
                % %                 faseApp=phase;
                % %
                % %                 X=[ sum(Amp.*xxx.^2) sum(Amp.*xxx.^1);
                % %                     sum(Amp.*xxx.^1) sum(Amp.*xxx.^0)];
                % %
                % %                 Y=[ sum(Amp.*faseApp.*xxx.^1) sum(Amp.*faseApp.*xxx.^0)]';
                % %
                % %                 Coeff=pinv(X)*Y;
                % %
                % %                 faseSim=Coeff(1).*xxx+Coeff(2); %considera solo la fase dello spettro
                % %                 %fine interpolazione mediante least mean squares
                % %
                % %             else %if realPart==0
                % %                 faseSim=phase; %considera la fase dello spettro+rumore
                % %             end

                % double Trapezoidal filter
                linear=1/(tt1(end)-tt1(1))*(tt1-tt1(1));
                Trap=ones(size(C));
                Trap(1:length(tt1))=linear;

                Trap(end:-1:end-length(tt1)+1)=linear; % Trapezoidal apodization

                yf_ok=(Dt.*(C.*Trap))*exps;
                SpettroC1(hh,:)=2*(yf_ok); %**Benedetto** 09/10/2020
                %%%**Benedetto** ha messo come commento le seguenti righe il 09/10/2020
                % %             SpettroC1(hh,:)=2*(yf_ok).*exp(-1j*faseSim);

            else %the OTHER CASES of apodization are treated here

                if center_find % **Benedetto** 14/01/2021
                    I=round([1:length(C)]*(C.^2)'./(sum(C.^2)+1)); % Strategy to find the baricenter (+1 is required to avoid dividing by 0)
                else
                    I=length(C)./2;
                end

                if trapz==0 %if there's no apodization
                    Trap=ones(1,length(C));
                    s_0=C; %row
                else %if there's apodization different from the trapezoidal one
                    Trap=Apodization(trapz,length(C),I,7);
                    s_0=C.*Trap; % row
                end

                % fourier computed with the exponential matrix and no loop involved
                yf_0=(Dt.*s_0)*exps;

                %%%**Benedetto** ha messo come commento le seguenti righe il 09/10/2020
                % %         phase=angle(yf_0); % This phase will be used to correct the spectral phase of the entire trace
                % %
                % %         if realPart %**Benedetto** 26/09/2020
                % %             %Interpolazione mediante least mean squares
                % %             xxx=f;
                % %             Amp=abs(yf_0).^2;
                % %             faseApp=unwrap(phase); % era senza unwrap
                % %
                % %             X=[ sum(Amp.*xxx.^2) sum(Amp.*xxx.^1);
                % %                 sum(Amp.*xxx.^1) sum(Amp.*xxx.^0)];
                % %
                % %             Y=[ sum(Amp.*faseApp.*xxx.^1) sum(Amp.*faseApp.*xxx.^0)]';
                % %
                % %             Coeff=pinv(X)*Y;
                % %
                % %             faseSim=Coeff(1).*xxx+Coeff(2); %considera solo la fase dello spettro
                % %             %fine interpolazione mediante least mean squares
                % %
                % %         else %if realPart==0
                % %                 faseSim=phase; %considera la fase dello spettro+rumore
                % %         end

                % %         SpettroC1(hh,:)=2*(yf_0).*exp(-1j*faseSim);

                SpettroC1(hh,:)=2*(yf_0); %**Benedetto** 09/10/2020


                % figure(30+hh);
                figure(h0); axes(h0b);
                plot(tt,C0+Averages(y,x),'m',tt,C,'r',tt,C.*Trap,'b',tt,Trap*max(C.*Trap),'k','linewidth',2); axis tight

                figure(h0); axes(h0c);
                plot(f,abs(SpettroC1(hh,:)),'b','linewidth',2); %real-->abs **Benedetto** 09/10/2020

                pause; [x,y,button]=ginput(1);

            end;

            sp=hh;

            % end;

            for hh=1:sp
                SpettroC2(hh,:)=SpettroC1(hh,:)./max(SpettroC1(hh,:)); % era abs(reference);
            end;

            % Frequency calibration
            h3=figure(3);

            subplot(2,1,1)
            plot(f,abs(SpettroC1),'linewidth',2); axis tight; %real-->abs **Benedetto** 09/10/2020
            title('Not normalized')
            legend;


            subplot(2,1,2)
            plot(f,abs(SpettroC2),'linewidth',2); axis tight; %real-->abs **Benedetto** 09/10/2020
            legend;

            h3a=gca;
            title('Normalized to white - Select Peaks for frequency calibration - RIGHT click when finished');

            figure(h3); axes(h3a);
            cal=0;
            [X,Y,BUTTON] = ginput(1);
            while BUTTON<=1
                cal=cal+1;
                answer=inputdlg('Corresponding wavelength (nm)','Wavelength calibration',1,{'500'});
                fr_realC(cal)=c./(str2double(answer)*1e-9)/1e12; % Frequency in THz
                fr_pseudoC(cal)=X;

                figure(h3); axes(h3a);
                [X,Y,BUTTON] = ginput(1);
            end

            file_totCal=[];

            % To be moved
            switch cal

                case {0,1}
                    load_cal=round(Dinput('\n\n    Load calibration [0: no - 1: yes]? ',0));

                    if load_cal

                        [filename_Cal, pathname_Cal] = uigetfile('*.mat', 'Load Calibration',dir0);
                        file_totCal=[pathname_Cal,filename_Cal];

                        stringa=[' *** Loaded calibration ',file_totCal,' ***'];
                        fprintf('\n\n%s',stringa);

                        load(file_totCal); % Contains  fr_real0 and f_0
                        cal=2; % For next steps, this means that calibration is available

                        fr_real=interp1(f_0,fr_real0,f);

                        figure(h3);

                        subplot(2,1,1);
                        plot(fr_real,abs(SpettroC1),'linewidth',2); axis tight %real-->abs **Benedetto** 09/10/2020
                        xlabel('Frequency (THz)')

                        % Correction frequency -> wavelength: SpettroC1*fr_real.^2/c
                        subplot(2,1,2);
                        plot(c./(fr_real*1e12)/1e-9,...
                            abs(SpettroC1.*fr_real(ones(1,size(SpettroC1,1)),:).^2/c),...
                            'linewidth',2); axis tight %real-->abs **Benedetto** 09/10/2020
                        xlabel('Wavelength (nm)')
                        % spectrumLabel(gca)

                    end
                case 2
                    M=(fr_realC(2)-fr_realC(1))/(fr_pseudoC(2)-fr_pseudoC(1));
                    Q=fr_realC(1)-M*fr_pseudoC(1);
                    fr_real=M*f+Q;

                    % Saving calibration
                    fprintf('\n  -- Saving calibration: calibration ranges --');

                    fr_real02=c/(abs(Dinput('\n\n    Shorter calibration wavelength (nm): ',300))*1e-9)/1e12; % highest frequency, in THz
                    fr_real01=c/(abs(Dinput('    Longer calibration wavelength (nm): ',1200))*1e-9)/1e12; % lowest frequency, in THz

                    f_01=(fr_real01-Q)/M;
                    f_02=(fr_real02-Q)/M;

                    f_0=linspace(f_01,f_02,500);
                    fr_real0=M*f0+Q;

                    figure(h3);

                    subplot(2,2,1);
                    plot(fr_real,abs(SpettroC1),'linewidth',2); axis tight %real-->abs **Benedetto** 09/10/2020
                    xlabel('Frequency (THz)')

                    % Correction frequency -> wavelength: SpettroC1*fr_real.^2/c
                    subplot(2,2,2);
                    plot(c./(fr_real*1e12)/1e-9,...
                        abs(SpettroC1.*fr_real(ones(1,size(SpettroC1,1)),:).^2/c),...
                        'linewidth',2); axis tight %real-->abs **Benedetto** 09/10/2020
                    % spectrumLabel(gca);
                    xlabel('Wavelength (nm)');

                    subplot(2,2,3);
                    plot(f,fr_real,f_0,fr_real0,fr_pseudoC,fr_realC,'*','linewidth',3);
                    xlabel('Pseudofrequency');
                    ylabel('Optical frequency [THz]');
                    legend('Interpolation1','Interpolation2','Measurement');

                    [filename_Cal, pathname_Cal] = uiputfile('*.mat', 'Save Calibration as',dir0);
                    file_tot=[pathname_Cal,filename_Cal];
                    save(file_tot,'f_0','fr_real0'); % f_0: pseudofrequency; f_real0: optical frequency (THz)

                    stringa=['  -- Calibration saved in ',file_tot,' --'];
                    fprintf('\n%s',stringa);

                otherwise
                    % polynomial interpolation
                    [P,S,MU] = polyfit(fr_pseudoC,fr_realC,2);
                    fr_real = polyval(P,f,[],MU);

                    fprintf('\n  -- Saving calibration: calibration ranges --');

                    fr_real02=c/(abs(Dinput('\n\n    Shorter calibration wavelength (nm): ',300))*1e-9)/1e12; % highest frequency, in THz
                    fr_real01=c/(abs(Dinput('    Longer calibration wavelength (nm): ',1200))*1e-9)/1e12; % lowest frequency, in THz

                    fr_real0=linspace(fr_real01,fr_real02,500);

                    [P1,S1,MU] = polyfit(fr_realC,fr_pseudoC,2);
                    f_0 = polyval(P1,fr_real0,[],MU);

                    figure(h3);

                    subplot(2,2,1);
                    plot(fr_real,abs(SpettroC1),'linewidth',2); axis tight %real-->abs **Benedetto** 09/10/2020
                    xlabel('Frequency [THz]')

                    % Correction frequency -> wavelength: SpettroC1*fr_real.^2/c
                    subplot(2,2,2);
                    plot(c./(fr_real*1e12)/1e-9,...
                        abs(SpettroC1.*fr_real(ones(1,size(SpettroC1,1)),:).^2/c),...
                        'linewidth',2); axis tight %real-->abs **Benedetto** 09/10/2020
                    xlabel('Wavelength [nm]')
                    % spectrumLabel(gca)

                    subplot(2,2,3);
                    plot(f,fr_real,f_0,fr_real0,fr_pseudoC,fr_realC,'*','linewidth',3);
                    xlabel('Pseudofrequency');
                    ylabel('Optical frequency [THz]');
                    legend('Interpolation1','Interpolation2','Measurement');

                    % Saving calibration
                    [filename_Cal, pathname_Cal] = uiputfile('*.mat', 'Save Calibration as',dir0);
                    file_totCal=[pathname_Cal,filename_Cal];
                    save(file_totCal,'f_0','fr_real0'); % f_0: pseudofrequency; f_real0: optical frequency (THz)

                    stringa=['  -- Calibration saved in ',file_tot,' --'];
                    fprintf('\n%s',stringa);
            end

            Dt=diff(tt);
            Dt(end+1)=Dt(end);

            % hhh=figure;
            map1=colormap(jet(512));
            % xx=linspace(0,255,256); %**Benedetto**
            % yy=linspace(0,63,512);
            % map1=interp1(xx,map,yy);
            %
            % zz=log(yy/63)-log(yy(2)/63); zz=zz/max(zz)*63; zz(1)=0;
            % maplog=interp1(xx,map,zz);
            % close(hhh);

        end
end

fprintf('\n  -- Image processing --');

Map_zero=Dinput('\n\n    Map Zeros and contrast (1 = yes, 0 = no): ',0);

num_fig=5;

if Map_zero

    conversionWindow=WConversion_cryst; %finestra di aiuto per le conversioni **BENEDETTO**
    set(conversionWindow,'OuterPosition', [205 -7 114.5714 35.1765]);

    fprintf('\n  NB: for better zero mapping, a narrow spectral band should be selected');

    fmin_zeros=Dinput('\n    Minimum pseudo-frequency: ',fmin);
    fmax_zeros=Dinput('    Maximum pseudo-frequency: ',fmax);
    resol_zeros=round(Dinput('    Number of spectral points: ',resol));

    close(conversionWindow);

    f_zeros=linspace(fmin_zeros,fmax_zeros,resol_zeros);
    exps_zeros=exp(-1i*2*pi*tt'*f_zeros);

    num_fig=num_fig+1;
    HyperMatrixFilt1=zeros(size(A));
    h=waitbar(0,'Finding Zeros and contrast');

    Zeros=zeros(a,b);
    Contrast=zeros(a,b);
    Peak=zeros(a,b);

    for y=1:a
        waitbar(y/a,h);

        for x=1:b

            C=squeeze(HyperMatrix(y,x,:))';

            % Delay: deduced from slope of phase
            Phase=unwrap(angle((Dt.*(C.*Trap))*exps_zeros));
            % P = polyfit(f,Phase,1);
            % Zeros(y,x)=-P(1)/(2*pi);

            Zeros(y,x)=-mean(diff(Phase)./diff(f_zeros))/(2*pi);
            Peak(y,x)=max(C);
            Contrast(y,x)=max(C)/Averages(y,x);

        end;

    end;

    close(h);

    center=round(size(Zeros)./2);
    Zeros=Zeros-Zeros(center(1),center(2));

    %%%%%%%%%%%

    filt=Dinput('    Smoothing level (0 = no smoothing): ',0);
    if filt>0
        H = fspecial('gaussian',10,filt/2);
        Zeros1 = imfilter(Zeros,H,'replicate'); % Gaussian filter for smoothing
    end;

    Zeros1 = Zeros; % Gaussian filter for smoothing

    %%%%%%%%%%%%%%%%%

    figure(num_fig+5);
    pcolor(Zeros); axis equal; axis off; shading interp;
    colormap(map1);
    hold on;
    % contour(Zeros1(:,:),[-10:2.5:10],'k'); %remember to comment this command if I don't want to see the contour lines
    hold off
    title('Zeros (steps)');
    colorbar;
    file_imageRGB=[file_tot(1:end-4),'_Zero_',num2str(num_fig+5)];
    print(file_imageRGB,'-djpeg','-r300');

    figure(num_fig+6);
    pcolor(Contrast); axis equal; axis off; shading interp;
    colormap(map1);
    title('Contrast');
    colorbar;
    Contrast_mean=mean(Contrast(:));
    Contrast_std=sqrt(var(Contrast(:)));

    fprintf(['\n  -- Contrast: ',num2str(Contrast_mean*100),' +/- ',...
        num2str(Contrast_std*100),' %% --\n\n']);

    file_imageRGB=[file_tot(1:end-4),'_Contrast_',num2str(num_fig+6)];
    print(file_imageRGB,'-djpeg','-r300');

    figure(num_fig+7);
    pcolor(Peak); axis equal; axis off; shading interp;
    colormap(map1);
    title('Interferogram amplitude');
    colorbar;
    file_imageRGB=[file_tot(1:end-4),'_Amplitude_',num2str(num_fig+6)];
    print(file_imageRGB,'-djpeg','-r300');

    Contrast1=Contrast(:);
    fprintf(['\n\n *** Contrast: ',...
        num2str(mean(Contrast1)*100),' +/- ',num2str(sqrt(var(Contrast1))*100),'%% ***']);

end;

%%%%%%%%%
% Eventuale filtro gaussiano:
%
% H = fspecial('gaussian',25,2);
% Zeros_out = imfilter(Zeros,H,'replicate');
%%%%%%%%%


Sp_filt=Dinput('\n\n    Generation of spectral Hypercube? (1 = yes, 0 = no): ',0);
num_fig=0;

if Sp_filt

    % parpool;
    num_fig=num_fig+1;

    % filt=Dinput('    Smoothing level (0 = no smoothing): ',0);

    fprintf('\n    Interferogram background subtraction:');
    linearFit=Dinput('\n          0 = constant line (faster) - 1 = linear trend (slower)- 2 =quadratic trend ? ',0);
    trapz=round(Dinput('\n    Apodization (0: none - 1: Happ-Genzel - 2: 3-term Blackman-Harris - 3: 4-term Blackman-Harris - 4: Triangular - 5: Supergaussian - 6: Trapezoidal)? ',1));
    center_find=round(Dinput('\n    Interferograms center calculation (0: exactly in the middle   - 1: function baricenter)? ',1));

    %%% **Martina** changed these dinput functions default values. They are the same as the initial ones, no matter of the eventual chosen
    %%% values in the preview (7.04.2025)

    if trapz==6 %**Benedetto**

        auto_max=Dinput('\n    Trapezoidal apodization. Retrieval of maximum: [0: automatic - 1: manual]',0);

        if auto_max==1
            hh=figure;
            plot(Corr_summ,'linewidth',2); axis tight;
            zoom xon; pause;
            [Center,~]=ginput(1);
            close(hh);

            I=round(Center);
            [maxC,~] = max(C);

            I1=find(tt-tt(I)<=(tt(I)-tt(1))); % tt(I)-tt(1) is the distance between the first
            % point and the peak. It will be the distance
            % between the peak and the last point
            if length(I1)<3, I1=[1 2 3]; end;

            tt1=tt(I1); % row

        end;

    end;

    %     if filt>0
    %         H = fspecial('gaussian',10,filt/2);
    %     end;

    HyperMatrixFilt1=zeros(size(A));

    %%% **Martina**
    if ~preview %calibration must be loaded now

        [filename_Cal, pathname_Cal] = uigetfile('*.mat', 'Load Calibration',dir0);
        file_totCal=[pathname_Cal,filename_Cal];

        stringa=[' *** Loaded calibration ',file_totCal,' ***'];
        fprintf('\n\n%s',stringa);

        load(file_totCal); % Contains  fr_real0 and f_0
        cal=2; % For next steps, this means that calibration is available

        Dt=diff(tt);
        Dt(end+1)=Dt(end);
    end



    if cal>1

        load 'RGB_Transmission' % Contains wl(nm), R, G, B transmission (range 0-1)

        conversionWindow=WConversion_cryst; %finestra di aiuto per le conversioni **BENEDETTO**
        set(conversionWindow,'OuterPosition', [205 -7 114.5714 35.1765]);

        f_low=Dinput('\n    Lower frequency (THz): ',0);
        f_high=Dinput('    Higher frequency (THz): ',c./min(wl*1e-9)/1e12);
        f_points=Dinput('    Number of spectral points: ',50);

        close(conversionWindow);

        % From real frequencies to pseudofrequencies
        f_low=interp1(fr_real0,f_0,f_low);
        f_high=interp1(fr_real0,f_0,f_high);

        f_1=linspace(f_low,f_high,f_points);
        fr_real=interp1(f_0,fr_real0,f_1);

        % Calculate RGB spectra for deducing the RGB color of the figure
        lambda_nm=c./(fr_real*1e12)/1e-9;
        [lambdaCIE, X_CIE, Y_CIE, Z_CIE] = colorMatchFcn('CIE_1931');
        X_nm=interp1(lambdaCIE, X_CIE,lambda_nm);
        Y_nm=interp1(lambdaCIE, Y_CIE,lambda_nm);
        Z_nm=interp1(lambdaCIE, Z_CIE,lambda_nm);

        R_nm=interp1(wl, R,lambda_nm);
        G_nm=interp1(wl, G,lambda_nm);
        B_nm=interp1(wl, B,lambda_nm);

    else
        conversionWindow=WConversion_cryst; %finestra di aiuto per le conversioni **BENEDETTO**
        set(conversionWindow,'OuterPosition', [205 -7 114.5714 35.1765]);
        f_low=Dinput('\n Lower frequency (pseudo): ',0);
        f_high=Dinput(' Higher frequency (pseudo): ',max(f));
        f_points=Dinput(' number of spectral points: ',50);
        close(conversionWindow);

        X_nm=ones(1,f_points);
        Y_nm=ones(1,f_points);
        Z_nm=ones(1,f_points);

        R_nm=ones(1,f_points);
        G_nm=ones(1,f_points);
        B_nm=ones(1,f_points);

    end;

    f=linspace(f_low,f_high,f_points);
    Df=(f_high-f_low)/(f_points-1);
    Deltaf=(f_high-f_low);

    exps=exp(-1i*2*pi*tt'*f);

    h=waitbar(0,'Generation of Hypercube');

    % Infos(a,b)=struct;
    Immagine_xyz=zeros(a,b,3);
    ImmagineRGB=zeros(a,b);
    Hyperspectrum_cube=zeros(a,b,f_points);

    Amplitude1=zeros(a,b);

    tic
    for y=1:a
        waitbar(y/a,h);

        for x=1:b % it was parfor

            C=squeeze(HyperMatrix(y,x,:))';



            if linearFit  % Subtracts a line instead of a constant

                FitPar = polyfit(tt,C,1);

                FitLinear = polyval(FitPar,tt);
                C=C-FitLinear;

            elseif linearFit==2

                FitPar = polyfit(tt,C,2);

                FitQuadr = polyval(FitPar,tt);
                C=C-FitQuadr;
            end

            %**Benedetto**
            if trapz==6 %the CASE OF TRAPEZOIDAL apodization is treated independently **Benedetto**

                if auto_max==1

                    [maxC,~] = max(C);

                    %I, I1, tt1: calculated before
                else

                    %I, I1, tt1: calculated for any time
                    [maxC,I] = max(C); % finds the position of the max

                    I1=find(tt-tt(I)<=(tt(I)-tt(1))); % tt(I)-tt(1) is the distance between the first
                    % point and the peak. It will be the distance
                    % between the peak and the last point
                    if length(I1)<3, I1=[1 2 3]; end;

                    tt1=tt(I1); % row

                end;
                y1=C(I1); % I select only the symmetric correlation

                minC = min(C); % finds the min

                maxC=maxC+Averages(y,x);
                minC=minC+Averages(y,x);

                Dt_1=diff(tt1);
                Dt_1(end+1)=Dt_1(end);

                exps_1=exp(-1i*2*pi*tt1'*f);

                s_1=y1.*Apodization(4,length(y1),I); % row

                % fourier computed with the exponential matrix and no loop involved
                yf_1=(Dt_1.*s_1)*exps_1;

                % double Trapezoidal filter
                linear=1/(tt1(end)-tt1(1))*(tt1-tt1(1));
                Trap_l=ones(size(C));
                Trap_l(1:length(tt1))=linear;

                Trap_l(end:-1:end-length(tt1)+1)=linear; % Trapezoidal apodization


                %%%%%%%%%%%%%%%% Era in fondo
                yf_ok=(Dt.*(C.*Trap_l))*exps;
                Spectrum=2*yf_ok;

                %%%**Benedetto** ha messo come commento le seguenti righe il 09/10/2020
                % %                     Spectrum=2*real((yf_ok).*exp(-1j*faseSim_1));

                Hyperspectrum_cube(y,x,:)=Spectrum;
                % Hyperspectrum_matrix(l_matrix,:)=Spectrum;

                HyperMatrixFilt1(y,x)=sum(Spectrum).*Df; % Integral of filtered spectrum
                Normalization=Averages(y,x)./HyperMatrixFilt1(y,x);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fine era in fondo

            else %the OTHER CASES of apodization are treated here

                if center_find % **Benedetto** 14/01/2021
                    I=round([1:length(C)]*(C.^2)'./(sum(C.^2)+1)); % Strategy to find the baricenter (+1 is required to avoid dividing by 0)
                else
                    I=length(C)./2;
                end

                if trapz==0 %if there's no apodization
                    Trap_1=ones(1,length(C));
                    s_0=C; %row

                else %if there's apodization different from the trapezoidal one
                    Trap_1=Apodization(trapz,length(C),I,7);
                    s_0=C.*Trap_1; % row

                end

                % fourier computed with the exponential matrix and no loop involved
                yf_0=(Dt.*s_0)*exps;

                Spectrum=2*yf_0;
                %%%**Benedetto** ha messo come commento le seguenti righe il 09/10/2020
                % %                 Spectrum=2*real((yf_0).*exp(-1j*faseSim_1));

                Hyperspectrum_cube(y,x,:)=Spectrum;

                HyperMatrixFilt1(y,x)=sum(Spectrum).*Df; % Integral of filtered spectrum
                Normalization=Averages(y,x)./HyperMatrixFilt1(y,x);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fine era in fondo

            end


        end;

    end;
    toc

    close(h);

    if Dinput('\n Save: [1] Absolute value; [0] Complex value: ',1)
        Hyperspectrum_cube=abs(Hyperspectrum_cube);
    end;

    %**Benedetto 23/06/2023 new way to save spectra hypercube dataset .mat
    %or .h5 (double, single, uint16) or mj2 lossless
    if exist('fr_real')
        save_hypercube(Hyperspectrum_cube,f,saturationMap,fr_real,file_totCal,[],dir0,file_tot);
    else
        save_hypercube(file_tot,Hyperspectrum_cube,f,saturationMap,[],[],[],dir0,file_tot);
    end

end

return;   %%% END OF PROGRAM
