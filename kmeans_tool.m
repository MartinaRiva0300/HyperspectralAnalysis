function [CLUSTERS,C]=kmeans_tool(Hyperspectrum_cube,n_clusters)
%kmeans_tool receives as input a spectral hypercube Hyperspectrum_cube and
%performes the k-means analysis by searching for a number of clusters equal
%to n_clusters

%kmeans_tool returns CLUSTERS, which is a 2D matrix cointaining in each
%pixel an integer value indicating the belonging cluster, and C which are
%the centroid spectra (C is a matrix where the column are the centroid
%spectra)

a=size(Hyperspectrum_cube,1);
b=size(Hyperspectrum_cube,2);
n_observations=a.*b; %number of observations, number of spectra in the image
p_variables=size(Hyperspectrum_cube,3); %number of variables, number of frequencies

X=reshape(Hyperspectrum_cube,[n_observations,p_variables]);
%the reshape is done in a columnwise: X is a list of vectors of dimension
%p_variables, the order of these vectors is by scanning the columns of the
%Hypespectrum_cube (after an entire column has been scanned, the nearby column is
%scanned)

[IDX,C]=kmeans(X,n_clusters);
CLUSTERS=reshape(IDX,a,b,1); %CLUSTERS is a 2D map where each pixels has an integer value indicating the belonging cluster
C=C'; %in this way the centroids are column spectra
end