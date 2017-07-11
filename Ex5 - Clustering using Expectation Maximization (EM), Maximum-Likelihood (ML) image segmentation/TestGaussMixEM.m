% function TestGaussMixEM
% 
% Randomly generates vectors in a three-dimensional feature space
% Calls EM-Algorithm to learn a statistical model based on this data
    
function TestGaussMixEM

    % Number of dimensions
    n_dims = 3;

    % number of desired components (cluster)
    n_comp = 8;

    % standard deviation of the clusters
    stddev_max = 0.3;

    % size (number of points) of each cluster
    cluster_size = 150;

    %feature vector initialization
    trainVect = zeros(n_comp,cluster_size,n_dims);
    % Mean value (for each cluster one vector of three components)
    % initialization
    Mean = randn(n_comp,n_dims);
    % Standard deviation initialization
    Stddev = randn(n_comp,n_dims)*stddev_max;

    % Random sampling of training vectors for each component according to 
    % the mean vector
    for i=1:n_comp
        tmp = zeros(cluster_size,1)+1;
        trainVect(i,:,:) = tmp*Mean(i,:);
        trainVect(i,:,:) = squeeze(trainVect(i,:,:)) + (tmp*Stddev(i,:)).*randn(cluster_size,n_dims);
    end

    % reshaping of vectors
    trainVect = reshape(trainVect,[n_comp*cluster_size,n_dims]);

    % GaussMixModel mittels EM lernen...
    model = LearnGaussMixModel(trainVect, n_comp);

end