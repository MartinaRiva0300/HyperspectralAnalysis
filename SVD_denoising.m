function [HyperMatrix_new]=SVD_denoising(HyperMatrix,t,f,image_ratio,num_SVs)

global threshold n_SVs

%function which performs the SVD denoising of the Temporal Hypercube
%(HyperMatrix)
%this SVD denoising algorithm is described in
%https://doi.org/10.1002/jrs.4824 (see supplementary)
%In such reference it is applied to spectral hypercubes, for this function
%it is applied to temporal hypercubes

%INPUTS
% - HyperMatrix : temporal hypercube
% - t : delay axis (motor positions)
% - f : pseudofrequencies axis (it is necessary for the spectrum preview)
% - image_ratio : ratio of the edge of the central part with respect to the
%                 image size (this is used to calculate the Spatial Signal
%                 Ratio figure of merit) [e.g. image_ratio = 0.25 takes 1/4
%                 of the image as central part]
% - num_SVs : number of first singular values (SVs) considered by the algorithm

%OUTPUTS
% - HyperMatrix_new : new denoised temporal hypercube

n_SVs=num_SVs;

HyperMatrix=double(HyperMatrix); %HyperMatrix could have been saved in single

if max(max(max(HyperMatrix)))<=0 %if HyperMatrix has been flipped (cross TWINS polarizers)
    hyp_flip=1;
else
    hyp_flip=0;
end

if hyp_flip==1
    HyperMatrix=-HyperMatrix;
end

edge_ratio=sqrt(image_ratio); %ratio between the edge of the central part 
%and the border of the image for the Spatial Signal Ratio calculation

threshold=image_ratio.*2; %starting threshold for Spatial Signal Ration to be considered

hhh=waitbar(1/2,'Performing SVD, please wait...');

tic;

[a,b,p_variables]=size(HyperMatrix);
n_observations=a.*b;

X=reshape(HyperMatrix,[n_observations,p_variables])'; %X: rows spectral axis, column spatial content

[U,S,V]=svds(X,n_SVs); %SVD (only first n_SVs Singular Values)

spatial_signal_ratio=zeros(1,size(V,2)); %spatial-signal ratio vector

a_start_mappa=round((a-a.*edge_ratio)./2);
b_start_mappa=round((b-b.*edge_ratio)./2);
a_end_mappa=a_start_mappa+round(a.*edge_ratio);
b_end_mappa=b_start_mappa+round(b.*edge_ratio);

%V columns represent the spatial weights for each spectral eigenfunction (PC)
%contained in the rows of U
for index=1:size(V,2) %index slides through the SVs (each associated to a column of V)
    
    SV_V_comp=reshape(V(:,index),[a,b]); %consider the index-th column of V as spatial image
    S_V_comp_FT=fft2(SV_V_comp); %take the spatial FT
    mappa=fftshift(abs(S_V_comp_FT));
    %in the FT map (mappa) the signal is located only in the center (low
    %spatial frequencies), the noise could be located everywhere
    Signal=mappa(a_start_mappa:a_end_mappa,b_start_mappa:b_end_mappa); %the central part is consider as Signal 
    Noise=mappa;
    Noise(a_start_mappa:a_end_mappa,b_start_mappa:b_end_mappa)=0; %the outside part is considered as only noise
    
    %define the "figure of merit" Spatial-Signal Ratio as the ratio between
    %the integral of the intensity in the Signal and in the Noise part of
    %the FT image
    spatial_signal_ratio(index)=sum(Signal(:))./sum(Noise(:)); %spatial-signal ratio associated to the index-th SV
end

time=toc;

waitbar(1,hhh,'Performing SVD, please wait...');
close(hhh);

fprintf('\nTime elapsed for SVD decomposition: %.1f s\n',time);

hhh=waitbar(1/2,'Opening the GUI, please wait...');

if length(spatial_signal_ratio)~=n_SVs %sometimes may happen
    n_SVs=length(spatial_signal_ratio);
end

hh=SVD_Spatial_signal_ratio(U,S,V,a,b,p_variables,spatial_signal_ratio, ...
threshold,n_SVs,t,f);

waitbar(1,hhh,'Opening the GUI, please wait...');
close(hhh);

%SVD_Spatial_signal_ratio directly modify the global variable 'threshold'

waitfor(hh); %while the SVD_Spatial_signal_ratio GUI is open pause the execution

%now 'threshold' is equal to the last value selected in the
%SVD_Spatial_signal_ratio GUI

SVs=diag(S);

SVs_swap=zeros(size(SVs));
SVs_swap(1:end-(length(SVs)-n_SVs))=SVs(1:n_SVs); %keep only the first n_SVs SVs

SVs=SVs_swap;
clear SVs_swap

SVs(spatial_signal_ratio<threshold)=0; %put to 0 all the SVs associated to low spatial_signal_ratio

S_new=diag(SVs); %recalculate the S matrix

X_new=(U*S_new*V')'; %recalculate the denoised dataset

HyperMatrix_new=reshape(X_new,a,b,p_variables); %reshape properly the dataset

if hyp_flip==1
    HyperMatrix_new=-HyperMatrix_new;
end

end