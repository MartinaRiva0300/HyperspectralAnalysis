%This program allows to select an AOI in the hyperspectral image. All the
%pixels outside/inside (out_NaN=1 / out_NaN=0) the AOI are put to Not a Number. The this program allows
%you to save the new masked Spectral Hypercube

out_NaN=1; %if out_Nan is 1 all the pixels out from the AOI are set to NaN, if out_Nan is 0 all the pixels in the AOI are set to NaN

dir0=('C:\Users\HARDi\Desktop\Tesi');

[filename_Hyper, pathname_Hyper] = uigetfile('*.mat', 'Load spectral Hypercube',dir0);
file_tot=[pathname_Hyper,filename_Hyper];

dir2=pathname_Hyper;

hh=waitbar(1/2,'Loading the Spectral Hypercube. Please wait...');

input_data=load(file_tot);

Hyperspectrum_cube=input_data.Hyperspectrum_cube;
fr_real=input_data.fr_real;
f=input_data.f;

if isfield(input_data,'saturationMap')==1 %if saturationMap has been saved together with Hyperspectrum_cube
    saturationMap=input_data.saturationMap;
else
    saturationMap=ones(size(Hyperspectrum_cube,1),size(Hyperspectrum_cube,2));
end

if isfield(input_data,'Intens')==1 %if Intens has been saved together with Hyperspectrum_cube **Benedetto** 12/07/2021
    Intens=input_data.Intens;
else
    Intens=sum(abs(Hyperspectrum_cube),3);
end

waitbar(2/2,hh);
close(hh);

figure;
imagesc(Intens);
M=linspace(0,1,512)';
MAP=[M M M];
colormap(MAP);
title('Please select the AOI mask');

[ThisAOI] = roipoly();

if out_NaN==1
    for index=1:length(fr_real)
        Frame=Hyperspectrum_cube(:,:,index);
        Frame(ThisAOI==0)=NaN;
        Hyperspectrum_cube(:,:,index)=Frame;
    end
else
    for index=1:length(fr_real)
        Frame=Hyperspectrum_cube(:,:,index);
        Frame(ThisAOI==0)=NaN;
        Hyperspectrum_cube(:,:,index)=Frame;
    end
end

Intens=sum(abs(Hyperspectrum_cube),3); %the Intens map is recalculated for the masked Hyperspectrum_cube

[filename_Spectra, pathname_Spectra] = uiputfile('*.mat', 'Save Current Hypercube',dir2);
            file_tot=[pathname_Spectra,filename_Spectra];
            
            h=waitbar(0.5,'Saving current Hypercube, please wait...');
            
            stringa=['  -- Saving current Hypercube in ',file_tot,' --'];
            fprintf('\n%s\n',stringa);
            
            % Checks the size of the matrix, before saving
            HyProp = whos('Hyperspectrum_cube') ;
            Gigabytes = HyProp.bytes/2^30;
            
            if Gigabytes<1.8
                save(file_tot,'f','fr_real','saturationMap','Intens','Hyperspectrum_cube');
            else
                save(file_tot,'f','fr_real','saturationMap','Intens','Hyperspectrum_cube','-v7.3');
            end;
            
            close(h);
            close all;
            clear all;