function Hyp_to_RamApp_generalPurpose(n_pixel_smoothing,n_wavenumber)

%code that loads a spectral hypercube in .h5, applies a smoothing (gaussian
%filter 5x5 pixels and stdev equal to n_pixel_smoothing) and save it in
%absolute uint16 .mat format readable with RamApp (array +
%Hypercube dataset)

    n_pixel_smoothing=round(n_pixel_smoothing);

    if n_pixel_smoothing<0
        n_pixel_smoothing=abs(n_pixel_smoothing);
    end

    %codice di H5HypercubeRead
        dir0=('C:\Users\HARDi\Desktop\Tesi\Codice immagini');

        [file_name,path_name] = uigetfile('*.h5', 'Load spectral Hypercube',dir0);
        filename=[path_name '\' file_name];

        info=h5info(filename); %struct with info of Datasets in H5 file

        if isempty(cell2mat({info.Datasets}))
            file_totCal=[];
        else
            H5_datasets_list=cell2mat({info.Datasets.Name}); %string listing all datasets in the H5 file

            if contains(H5_datasets_list,'file_totCal')
                file_totCal=cell2mat(h5read(filename,'/file_totCal')); %if calibration string exists load it
            else
                file_totCal=[];
            end
        end

        info_datasets=h5info(filename,'/SpectralHypercube'); %struct with info of Datasets in SpectralHypercube dataset in H5 file
        H5_variables_list=cell2mat({info_datasets.Datasets.Name}); %string listing all datasets (variables) names inside the SpectralHypercube dataset in the H5 file

        Hyperspectrum_cube=h5read(filename,'/SpectralHypercube/Hyperspectrum_cube');
        f=h5read(filename,'/SpectralHypercube/f');
        saturationMap=h5read(filename,'/SpectralHypercube/saturationMap');

        if contains(H5_variables_list,'fr_real')
            fr_real=h5read(filename,'/SpectralHypercube/fr_real'); %if calibrated frequencies axis exists load it
        else
            fr_real=[];
        end

        if contains(H5_variables_list,'Intens')
            Intens=h5read(filename,'/SpectralHypercube/Intens'); %if Intens matrix preview exists
        else
            Intens=[];
        end

        if contains(H5_variables_list,'minimum')
            minimum=h5read(filename,'/SpectralHypercube/minimum'); %if minimum (to reconvert uint16 dataset) exists load it
        else
            minimum=[];
        end

        if contains(H5_variables_list,'maximum')
            maximum=h5read(filename,'/SpectralHypercube/maximum'); %if maximum (to reconvert uint16 dataset) exists load it
        else
            maximum=[];
        end

        if isa(Hyperspectrum_cube,'uint16')
            Hyperspectrum_cube=double(Hyperspectrum_cube)./(2.^16-1).*maximum+minimum; %double reconversion procedure from uint16
            fprintf('\n The loaded hypercube was saved in uint16 type\n\n'); %message to explicit the hypercube format
        elseif isa(Hyperspectrum_cube,'single')
            Hyperspectrum_cube=double(Hyperspectrum_cube); %reconversion to double
            fprintf('\n The loaded hypercube was saved in single type\n\n'); %message to explicit the hypercube format
        else
            fprintf('\n The loaded hypercube was saved in double type\n\n'); %message to explicit the hypercube format
        end

        if size(Hyperspectrum_cube,3)==2.*length(f) %if the spectral hypecube is complex the real and imag are stacked one over the other
            Hyperspectrum_cube_real=Hyperspectrum_cube(:,:,1:size(Hyperspectrum_cube,3)./2);
            Hyperspectrum_cube_imag=Hyperspectrum_cube(:,:,size(Hyperspectrum_cube,3)./2+1:end);
            clear Hyperspectrum_cube
            Hyperspectrum_cube=zeros(size(Hyperspectrum_cube_real));
            Hyperspectrum_cube=Hyperspectrum_cube_real+1i.*Hyperspectrum_cube_imag;
            clear Hyperspectrum_cube_real Hyperspectrum_cube_imag
        end
    %fine codice H5HypercubeRead
    
    if exist('n_wavenumber')==0
        n_wavenumber=size(Hyperspectrum_cube,3);
    end

    if n_pixel_smoothing>0
        %H = fspecial('average',n_pixel_smoothing)
        %'average' è un filtro rettangolare
        %era H = fspecial('gaussian',10,filt/2);
        %filt/2 era la st.dev. della gaussiana (filt
        %inserito dall'utente)
        H = fspecial('gaussian',5,n_pixel_smoothing./2);

        %il filtro non deve necessariamente essere
        %rettangolare, può essere anche Gaussiano o altro
        %(sono tutti buoni filtri)

        sp_points=size(Hyperspectrum_cube,3);

        h=waitbar(0,'Smoothing Hypercube, please wait...');

        for sp=1:sp_points

            waitbar(sp/sp_points,h);
            Hyperspectrum_cube(:,:,sp)=imfilter(Hyperspectrum_cube(:,:,sp),H,'replicate'); % imfilter smoothing (it works properly even with complex hypercubes)
        end

        close(h);
    end
    
    c=299792458; %speed of light [m/s]
    axisType=Dinput('Select the type of x axis:\n1) wavenumber; \n2)frequency (Hz); \n3) wavelength (nm)\n', 3);
    switch axisType
        case 1
            wavenumber=(c./532.*1e-3-fr_real)./c.*1e12./100;
            wavenumber=flip(wavenumber);

            Hyperspectrum_cube=abs(Hyperspectrum_cube);
            Hyperspectrum_cube=flip(Hyperspectrum_cube,3);

        case 2
            wavenumber=fr_real*10^(12); %in Hz

        case 3
            wavenumber=c./fr_real*10^(-3);
            wavenumber=flip(wavenumber);

            Hyperspectrum_cube=abs(Hyperspectrum_cube);
            Hyperspectrum_cube=flip(Hyperspectrum_cube,3);


            if n_wavenumber~=length(wavenumber)

                wavenumber1=linspace(min(wavenumber),max(wavenumber),n_wavenumber);

                a=size(Hyperspectrum_cube,1);
                b=size(Hyperspectrum_cube,2);
                Hyperspectrum_cube1=zeros(a,b,n_wavenumber);

                hhhh=waitbar(0,'Reinterpolation of the Hypercube, please wait...');
                for yy=1:a
                    for xx=1:b
                        spectrum=squeeze(Hyperspectrum_cube(yy,xx,:));
                        Hyperspectrum_cube1(yy,xx,:)=interp1(wavenumber,spectrum,wavenumber1);
                        waitbar((xx+(yy-1)*b)./(a.*b),hhhh);
                    end
                end
                close(hhhh);

                clear wavenumber spectrum
                wavenumber=wavenumber1;
                clear wavenumber1
                clear Hyperspectrum_cube
                Hyperspectrum_cube=Hyperspectrum_cube1;
                clear Hyperspectrum_cube1

            end

            maximum=max(max(max(real(Hyperspectrum_cube))));
            Hyperspectrum_cube=Hyperspectrum_cube./maximum;
            Hyperspectrum_cube=uint16((2.^16-1) * mat2gray(Hyperspectrum_cube,[0 1]));

            if not(isfolder([path_name '\RamApp_Data']))
                mkdir([path_name '\RamApp_Data']); %create RamApp_Data folder if it doesn't exist
            end

            file_new=[path_name '\RamApp_Data\' file_name(1:end-3) '_SMOOTH' num2str(n_pixel_smoothing) '_RamApp.mat'];

            h=waitbar(1/2,'Saving Hypercube, please wait...');
            save(file_new,'wavenumber','Hyperspectrum_cube');
            waitbar(1,h);
            close(h);




    end

