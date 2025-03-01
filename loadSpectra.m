function [wavenumber,spectra,spectra_Norm,wavelength,spectra_wavelength,frequency,spectra_frequency]=loadSpectra(indexesOfSpectra)

%loadSpectra
%input:  -indexesOfSpectra, the indexes of the spectra in the order of
%         selection in HyperspectralAnalysis_Spectrum; the spectra must be 
%         saved in the txt format decided in HyperspectralAnalysis_Spectrum
%
%output: -wavenumber, the vector of the wavenumbers in cm^(-1) calculated 
%         with respect to the laser wavelength
%        -spectra, the spectra that you want organized in columns
%        -spectraNorm, the spectra that you want normalized with respect to
%         the maximum of each one
%        -wavelength, the vector of the wavelength in nm
%        -spectra_wavelength, the spectra in frequency and Jacobian
%         corrected

c=299792458; %speed of light
dir0=('C:\Users\HARDi\Desktop\Tesi\Codice immagini');
[filename, pathname] = uigetfile('*.txt', 'Load saved spectra',dir0);
file_tot=[pathname,filename];

input_data=load(file_tot);

prompt={'Laser wavelength (nm)'};
name='Write laser wavelength';
numlines=1;
defaultanswer={'532'}; %era 531.4 per il Verdi

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer=inputdlg(prompt,name,numlines,defaultanswer,options);
laser=str2num(cell2mat(answer));

laser_frequency_THz=c/(laser*10^(-9))*10^(-12);

wavenumber=(-input_data(:,1)+laser_frequency_THz)./c.*10^(12)./100; %wavenumber in cm^(-1)
spectra=input_data(:,1+(indexesOfSpectra));
spectra_Norm=spectra./max(spectra);
%Jacobian correction
wavelength=c./input_data(:,1).*10^(-3);
spectra_wavelength=spectra.*wavelength.^(2)./c;
%spectra in frequency as they have been saved in the txt file
frequency=input_data(:,1);
spectra_frequency=spectra;

end