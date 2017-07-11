    
function ApplyGaussMixEM
    
    % read image with 3 channels!
    [file, path, image, img_size] = read_image('Select input file for segmentation');
 
    % Number of dimensions
    n_dims = img_size(3);
    
    %--------------------------------------------------------------------------
    % number of desired components (clusters)
    % vary this parameter to find an appropriate value for the input
    % image (Task )!!!
    n_comp = 8;
    %--------------------------------------------------------------------------
    
    % reshaping of vectors for input of EM
    trainVect = reshape(image, [img_size(1)*img_size(2),n_dims] );
    
    % sample the vectors to reduce amount of data
    desired_number = 1000;
    step = img_size(1)*img_size(2) / desired_number; 
    indices = round(1:step:img_size(1)*img_size(2));
    trainVect = trainVect(indices,:);
    
    % filter out edge points (lead to covariance-matrices which are not invertible)
    s_t = sum(trainVect,2);
    test = s_t ~= 3.0;
    trainVect = trainVect(test, :);
    
    % GaussMixModel mittels EM lernen...
    model = LearnGaussMixModel(trainVect, n_comp);

    % classify pixels
    ClassImg = ClassifyImage(model, image);

    % visualize result
    figure(2); 
    subplot(1,2,1), imshow(image), title('Original Image');
    subplot(1,2,2), imshow(ClassImg,[]), colormap(jet), title('Classification Result');
end



% PASTE HERE YOUR IMPLEMENTED FUNCTION CalcLnVectorProb (TASK B.a)
function LnVectorProb = CalcLnVectorProb(model, trainVect)

% IMPLEMENT THIS FUNCTION (TASK A.a)
temp = size(model.weight); % to evade the matlab problem of creating a temp array
Nc = temp(1);
temp2 = size(trainVect);
Nx = temp2(1);
LnVectorProb = zeros(Nc,Nx);

for c=1:Nc
    CoVMat = squeeze(model.covar(c,:,:));
    CovDet = det(CoVMat);
    CovInv = inv(CoVMat); %matlab suggest never to multiple by inverse
   % matrix instead to use CoVMat/Trainshift for linear equations
    for x=1:Nx
        TrainShift= trainVect(x,:) - model.mean(c,:);
        LnVectorProb(c,x) = log(model.weight(c)) - 0.5*(log(CovDet) + TrainShift*CovInv*TrainShift');
    end
end
end




%--------------------------------------------------------------------------
function  ClassImg = ClassifyImage(model, image)
    
    % image dimensions
    s = size(image);
   
    % reshaping of feature vectors
    FeatureVects = reshape(image, [s(1)*s(2),s(3)]);
    
    % probability of all vectors for all clusters (classes)
    LnVectorProb = CalcLnVectorProb(model, FeatureVects);  
    
    % get the maximum value --> this is the corresponding class membership
    [max_values, max_pos]  = max(LnVectorProb,[],1); 
    
    % reshape vector to result array
    ClassImg = uint8(reshape(max_pos, s(1:2)));
end

%--------------------------------------------------------------------------
function [file, path, image, s] = read_image(text)

    % open a dialogue to pick afile
    [file, path] = uigetfile('*.*', text);
 
    % read image and convert to double [0,...,1]
    image = mat2gray( imread([path,file]) );
    
    % determine size/dimensions of image
    s = size(image);
end