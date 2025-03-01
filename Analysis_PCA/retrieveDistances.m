function dist = retrieveDistances(method)
% This function creates a list of distances accepted by the clustering
% algorithm specified by the input method.
method = lower(method);
switch method
    case 'k-means'
        dist = {'sqeuclidean','cityblock','cosine','correlation','hamming' };
    case 'dbscan'
        dist = {'euclidean', 'squaredeuclidean', 'seuclidean', 'precomputed', ...
            'mahalanobis', 'cityblock','minkowski', 'chebychev', 'cosine', ...
            'correlation', 'hamming', 'jaccard', 'spearman', '@distfun'};
    case 'hierarchical'
        dist = {'euclidean', 'squaredeuclidean', 'seuclidean', ...
            'fasteuclidean', 'fastsquaredeuclidean', 'fastseuclidean', ...
            'mahalanobis', 'cityblock','minkowski', 'chebychev', 'cosine', ...
            'correlation', 'hamming', 'jaccard', 'spearman', '@distfun'};
    case 'clusterdata'
        dist = {'euclidean', 'squaredeuclidean', 'seuclidean', ...
            'mahalanobis', 'cityblock','minkowski', 'chebychev', 'cosine', ...
            'correlation', 'hamming', 'jaccard', 'spearman', '@distfun'};
    otherwise
        warndlg('Method not recognized')
        dist = {};
end