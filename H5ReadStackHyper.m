function [data,el_size,step] = H5ReadStackHyper(filename,ROI)
% el_size = [dx,dy,dz]
% step = [start_pos,step,step_num]
location = '/measurement/hyper/t0/c0/';
dataset = 'image';
h5disp(filename,location)
if nargin<2
    data = h5read(filename,[location,dataset]);% occhio a invertire X e Y dell'array ROI
else
    data = h5read(filename,[location,dataset],[ROI(2) ROI(1) 1],[ROI(4) ROI(3) inf]);
end
el_size = h5readatt(filename,[location,dataset],'element_size_um');
step(1) = h5readatt(filename,['/measurement/hyper/settings'],'start_pos');
step(2) = h5readatt(filename,['/measurement/hyper/settings'],'step');
step(3) = h5readatt(filename,['/measurement/hyper/settings'],'step_num');
end

