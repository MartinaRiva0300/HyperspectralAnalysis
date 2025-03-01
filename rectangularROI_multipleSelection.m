%% Rectangular ROI selection
%this script allows one to select more than one rectangular ROI with the
%same aspect ratio; it saves the selected areas mean spectra and standard
%deviations in a .txt file and the image with the selected areas in .jpeg

[filename_Hyper, pathname_Hyper] = uigetfile('*.mat;*.h5;*.mj2', 'Load spectral Hypercube'); 
file_tot=[pathname_Hyper,filename_Hyper];
input=load(file_tot);
f=input.f;

Hyperspectrum_cube=input.Hyperspectrum_cube;


fr_real=input.fr_real;
c=299792458;

a=size(Hyperspectrum_cube, 1);
b=size(Hyperspectrum_cube, 2);
if size(Hyperspectrum_cube,3)==2.*length(f) %if the spectral hypecube is complex the real and imag are stacked one over the other
    Hyperspectrum_cube_real=Hyperspectrum_cube(:,:,1:size(Hyperspectrum_cube,3)./2);
    Hyperspectrum_cube_imag=Hyperspectrum_cube(:,:,size(Hyperspectrum_cube,3)./2+1:end);
    clear Hyperspectrum_cube
    Hyperspectrum_cube=zeros(size(Hyperspectrum_cube_real));
    Hyperspectrum_cube=Hyperspectrum_cube_real+1i.*Hyperspectrum_cube_imag;
    clear Hyperspectrum_cube_real Hyperspectrum_cube_imag
end
fig1=figure;
imagesc(squeeze(abs(Hyperspectrum_cube(:,:,250)))); colormap gray;
roi=drawrectangle();
% extract top-left corner position, width and height
rectPosition=roi.Position;
% Extract x, y, width, height
%x = rectPosition(1);
%y = rectPosition(2);
width = rectPosition(3);
height = rectPosition(4);

% Calculate the aspect ratio
aspectRatio = width / height;
%display(aspectRatio)
%Keep the aspect ratio and select the new coordinates: clicking in a point
%without drawing any rectangle I can derive the coordinates of the top-left
%corner while I use the previous aspect ratio
mask=createMask(roi);
%mask3d=repmat(mask, [1,1, length(fr_real)]);
count=0;
for y=1:a
    for x=1:b
        if mask(y,x)
            count=count+1;
            spectrum(count, :)=squeeze(Hyperspectrum_cube(y,x,:));
        end
    end
end
spectrum_ave(1,:)=abs(mean(spectrum));
std_spectrum(1,:)=sqrt(var(abs(spectrum)));


newRoi=Dinput('\n\n   Select new ROI with same aspect ratio: (0 = no, 1 = yes): ',1);
i=2;
while newRoi
    [x, y] = ginput(1);
    roi=drawrectangle('Position', [x, y, width, height]);
    mask=createMask(roi);
    mask3d=repmat(mask, [1,1, length(fr_real)]);
    count=0;
    for y=1:a
        for x=1:b
            if mask(y,x)==1
                count=count+1;
                spectrum(count, :)=squeeze(Hyperspectrum_cube(y,x,:));
            end
        end
    end
    spectrum_ave(i,:)=abs(mean(spectrum));
    std_spectrum(i,:)=sqrt(var(abs(spectrum)));
    newRoi=Dinput('\n\n   Select new ROI with same aspect ratio: (0 = no, 1 = yes): ',1);
    i=i+1;
end

fig2=figure;
for j=1:(i-1)
    plot(c./fr_real, spectrum_ave(j,:), c./fr_real, spectrum_ave(j,:)+std_spectrum(j,:), c./fr_real, spectrum_ave(j,:)-std_spectrum(j,:))
    hold on
end


% Saving selected spectra
[filename_Spectra, pathname_Spectra] = uiputfile('*.txt', 'Save Spectra as',pathname_Hyper);
file_tot=[pathname_Spectra,filename_Spectra];

spectra_export=double([fr_real' spectrum_ave' std_spectrum']);

stringa=['  -- Saving Selected Spectra in ',file_tot,' --'];
fprintf('\n%s\n',stringa);
fprintf(['\n  Data format:\n\n    frequency[THz] (1 col)  Intensity[a.u.] Intensity_std[a.u.]\n\n']);
save(file_tot,'spectra_export','-ASCII'); % Selected Spectra

saveas(fig1,[file_tot(1:end-3),'jpg'],'jpeg');
