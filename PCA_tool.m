function [COEFF,SCORE]=PCA_tool(Hyperspectrum_cube)
%PCA_tool receives as input a spectral hypercube Hyperspectrum_cube and
%performes the PCA analysis

%PCA_tool returns COEFF whose columns are the PCs and SCORE which is the
%hypercube as a function of the PCs (the x and y coordinates are the x and
%y coordinates of the hyperspectral image, the z coordinate is the PC
%coordinates that now replace the frequencies)

a=size(Hyperspectrum_cube,1);
b=size(Hyperspectrum_cube,2);
n_observations=a.*b; %number of observations, number of spectra in the image
p_variables=size(Hyperspectrum_cube,3); %number of variables, number of frequencies

X=reshape(Hyperspectrum_cube,[n_observations,p_variables]);
%the reshape is done in a columnwise: X is a list of vectors of dimension
%p_variables, the order of these vectors is by scanning the columns of the
%Hypespectrum_cube (after an entire column has been scanned, the nearby column is
%scanned)

[COEFF,PCA_SCORE]=pca(X);
SCORE=reshape(PCA_SCORE,a,b,p_variables);

end