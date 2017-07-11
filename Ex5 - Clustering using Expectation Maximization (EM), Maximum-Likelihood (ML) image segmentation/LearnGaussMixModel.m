
% Implementation of a split-based EM-algorithm
%---------------------------------------------
%
% Inputs:
%   trainVect: array of size n_featurevects x n_dimensions
%   n_comp: number of desired components in the model
%
% Outputs: 
%   model: structure with estimated model
%   model.weight: (n_components x 1) vector with weight for each component
%   model.mean: (n_components x n_dims=3) matrix with mean vectors for each
%                component --> mean of vector i: model.mean(i,:)
%   model.covar: (n_components x n_dims=3 x n_dims=3) matrix with
%                covariance matrices for each component
%                squeeze(model.covar(i,:,:)) returns the i'th covar matrix
%
% main function of the EM-algorithm
function model = LearnGaussMixModel(trainVect, n_comp)

    % initialization of the model using a structure
    % at the starting point the algorithm will always be initialized with 
    % one cluster 
    model.weight(1,:)=1;
    model.mean=[0,0,0];
    model.covar(1,:,:)=[1 0 0; 0 1 0; 0 0 1];
    
    % threshold for stopping the iteration
    eps=10^-6;
    
    % loop over the desired number of components
    for i=1:n_comp
        
        % the first overall model probability is -infinity
        LastPX=-inf;
        
        % calculate the logarithmic overall model probability
        LnTotalProb = CalcLnTotalProb(model, trainVect);
        
        % while the threshold is bigger than the difference of the overall
        % probabilities from this and the last iteration...
        while (LnTotalProb-LastPX>eps)
            LastPX = LnTotalProb;
            
            % E-step:
            % compute for each feature vector the probabilities for each
            % component
            LnCompProb = GmmEStep(model, trainVect);
            
            % M-step:
            % Maximize the model by reestimating the model parameters using
            % the probabilities of E-step
            model = GmmMStep(model, trainVect, LnCompProb);
            
            % again compute the overall probability of the model
            % since EM always converges, this value is always higher than
            % in the last iteration
            LnTotalProb = CalcLnTotalProb(model, trainVect);      
        end
        
        % clear current figure window
        clf
        % plot the estimated GMM
        PlotGMM(model,trainVect);
        % flush
        drawnow;
        
        % find a component to split into two and init them
        % but only, if the desired number n_comp is not reached
        if i < n_comp
            model = InitNewComponent(model, trainVect);
        end
    end
end


%--------------------------------------------------------------------------
% logarithmic probability of all vectors for all components
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
% E-Step:
% calculation of the probabilities of each feature vector wrt all existing
% components using the current model parameters
function LnCompProb = GmmEStep(model, trainVect)
% IMPLEMENT THIS FUNCTION (TASK A.b)
Nc = numel(model.weight);
LnVectorProb = CalcLnVectorProb(model,trainVect);
MaxLVP = max(LnVectorProb, [],1);
ScalFactor= repmat(MaxLVP, Nc,1);
Denom = MaxLVP +log(sum(exp(LnVectorProb - ScalFactor)));
Denom = repmat(Denom,Nc,1);
LnCompProb = LnVectorProb - Denom;
end

%--------------------------------------------------------------------------
% M-Step:
% Estimation of new model parameters according the calculated probabilities
% of the E-Step
function model = GmmMStep(model, trainVect, LnCompProb)
% IMPLEMENT THIS FUNCTION (TASK A.c)   
    LCPExp = exp(LnCompProb);
    [Nv2,~] = size(trainVect);
    Ni = length (model.weight);
    Np = zeros(Ni,Nv2);
    SigmaCa = zeros(Ni,Nv2);
    for z=1:Ni
        for u=1:Nv2
            TrainShift = trainVect(u,:) - model.mean(z,:);
            Np = LCPExp(z,u);
            SigmaCa = Ni*Np*(TrainShift*TrainShift');
        end
    end
    AlphaC = Np/Nv2;
    SigmaC = SigmaCa/Np;
%     model = SigmaC;
    
    



end

%--------------------------------------------------------------------------
% calculation of the global probability given the current model and the
% feature vectors
function LnTotalProb = CalcLnTotalProb(model, trainVect)
    
    % get the current number of components in the model
    n_comp = numel(model.weight);

    % logarithmic probability for all vectors in all components
    LnVectorProb = CalcLnVectorProb(model, trainVect);
        
%     % the log of a sum cannot easity be computed from single log-values!
%     % so for this step we have to use the exp-function and afterwards take
%     % the log of the sum! (log of a product-->sum log values, but there is 
%     % no such rule for log of a sum!)
%     s = sum(exp(LnVectorProb),1);
%     LnTotalProb = sum(log(s));

    % the result abuve could be wrong tue to very small values after exp...
    % to be safe, we can compute LnTotalProb using a scaling factor:
    %---------------------------------------------------------------------
    % use scaling factor c = max of log values
    % wrt a feature vector --> one scale factor for each feature vctor

    % max probability for each feature vector
    % [ the maximum probability tells us, to which cluster each vector
    % belongs]
    % this value is used for scaling the probabilities in order to avoid
    % numerical problems for computation of the sum
    max_LnVectorProb = max( LnVectorProb,[],1 );
    
    % resize this array to size of LnTotalProb
    scaling_factors = repmat(max_LnVectorProb, n_comp, 1);
    
    % scaling of logarithmic probabilities before using exp in order to
    %  avoid numerical problems:
    % 1) subtract scaling_factors from LnVectorProb (scaling)
    % 2) take exp of the result (should be no problem after scaling)
    % 3) summarize the n_comp values for each feature vector (as desired)
    % 4) take the logarithm of the sums
    % 5) add the maximum to the result of 4 ("unscaling")
    LnVectorProb_new = max_LnVectorProb + log(sum(exp(LnVectorProb - scaling_factors),1));
 
    % sum all log values to get global model probability
    LnTotalProb = sum(LnVectorProb_new);
end




%--------------------------------------------------------------------------
% function for plotting a result of EM-estimation
% 
% Inputs:
%   model: Gaussian Mixture Model Parametes (structure)
%   features: feature vectors
% 
%  - plots a mixture model in a 3dim- plot (feature vector has to be 3-dimensional!)
%  - Feature vectors are plotted as green dots
%  - Means of components: red circles
%  - Covariance matrices: plotted using the three main axes of the
%    ellipsoid
function PlotGMM(model, trainVect)

    % number of components in the model
    n_comp = numel(model.weight);
    n_dims = size(trainVect, 2);
    
    hold on;
    % plot feature vector points
    plot3(trainVect(:,1),trainVect(:,2),trainVect(:,3), 'g.','MarkerSize',7);

    % plot elements of the estimated components:
    for i=1:n_comp
        
        % eigenvektor / eigenwert decomposition
        [eVec,eVal] = eig(squeeze(model.covar(i,:,:)));

        % plotting of mean values
        mean = squeeze(model.mean(i,:));
        plot3(mean(1),mean(2),mean(3),'ro');

        % derivation and plotting of the three main axes of the cvariance
        % matrices
        for i=1:n_dims
            devVec = (sqrt(eVal(i,i)) * eVec(:,i))*[-1,1];
            plot3(mean(1) + devVec(1,:), mean(2) + devVec(2,:), mean(3) + devVec(3,:),'b');
        end
    end

    % rotate 3D view and setting of title
    hold off;
    view([19,25]);
    grid('on');
    title(['Gaussian Mixture Model (',num2str(n_comp),' components)']);
end

% ------------------------------------------------------------------------
% function for splitting a single component and initialization 2 new ones
%
% Inputs:
%
%   model: GMM parameters
%   features: feature vectors
%   features(j,:): feature vector of a single pixel
%   newModel: updated model with new initialized component
%
% adds a new component (cluster) to the current model
%
% --> analyzes the current model and identifies the weakest component
% --> weak: doesn't fit well to the corresponding feature vectors
% the weakest component will be splitted into two new ones
function NewModel = InitNewComponent(model, trainVect)

    % number of components
    n_comp = numel(model.weight);

    % Number of dimensions (shall be three here!)
    n_dims = size(trainVect, 2);

    % the biggest component will be splitted to get a balanced size of
    % components --> not the optimal criterium!!!!
    % size corresponds to weights...
    [ignore, splitComp] = max(model.weight);

    % calculate new weight vector, mean and covariance
    newWeight = zeros(n_comp+1,1);
    newMean = zeros(n_comp+1,n_dims);
    newCovar = zeros(n_comp+1,n_dims,n_dims);

    % copy old values into new arrays
    newWeight(1:n_comp) = model.weight;
    newMean(1:n_comp,:) = model.mean;
    newCovar(1:n_comp,:,:) = model.covar;

    % Component splitComp will be splitted along its dominant axis
    [eVec,eVal] = eig(squeeze(newCovar(splitComp,:,:)));
    [ignore, majAxis] = max(diag(eVal));
    devVec = sqrt(eVal(majAxis,majAxis)) * eVec(:,majAxis)';

    % initialize new component
    % --> half of the points --> half weight
    newWeight(n_comp+1) = 0.5*newWeight(splitComp);
    % shift new mean to half of length along dominant axis
    newMean(n_comp+1,:) = newMean(splitComp,:) - 0.5*devVec;
    % make covariance a little bit smaller
    newCovar(n_comp+1,:,:) = newCovar(splitComp,:,:) / 4.0;

    % update also the (old) splitted component
    % also half of the points
    newWeight(splitComp) = newWeight(n_comp+1);
    % shift comonent center to other direction along dominant axis
    newMean(splitComp,:) = newMean(n_comp+1,:) + devVec;
    % take same smaller covariance matrix
    newCovar(splitComp,:,:) = newCovar(n_comp+1,:,:);

    % store new parameters in model
    NewModel.weight = newWeight;
    NewModel.mean = newMean;
    NewModel.covar = newCovar;
end
