%This program converts the .h5 file into a .mat file that can be feeded to
%HypespectralAnalysis_Time.m. The .mat temporal hypercube will be saved in
%the same folder of the .h5 file. The .mat temporal hypercube is a struct
%file with HyperMatrix (the temporal hypercube) and t (the delay axis)

%dir0=('C:\Users\HARDi\OneDrive\Desktop\Tesi');

[filename_Hyper, pathname_Hyper] = uigetfile('*.h5', 'Load Hypercube',dir0);

path(path,pathname_Hyper);

[HyperMatrix,el_size,step]=H5ReadStackHyper(filename_Hyper);
% el_size = [dx,dy,dz]
% step = [start_pos,step,step_num]
%start_pos in mm
%step in um

t=step(1).*1e3:step(2):step(1).*1e3+step(2).*(step(3)-1); %t is in um

HyperMatrix=double(HyperMatrix); %HyperMatrix is converted from uint16 to double
HyperMatrix=permute(HyperMatrix,[2 1 3]); %change X and Y as they are in the measurement setup

filename=[filename_Hyper(1:end-3) '_hyp.mat'];
file_tot=[pathname_Hyper,filename];

h=waitbar(0.5,'Saving hypercube in .mat file');
save(file_tot,'t','HyperMatrix','-v7.3');
close(h);

clear all;