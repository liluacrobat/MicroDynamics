function [ACC, w_history, Predict_value] = script_crossvail_LOGO(X, label, Para)
%% ==========================================================================
% Evaluate the performance of LOGO in k-fold cross-validation
%
%--------------------------------------------------------------------------
% INPUT:
%     X     :  Training data: [x1,x2,...xn] Each column is an observation
%     label :  Class label = {1,2,...,C}
%     Para  :  Parameters 
%         -- distance
%              Distance (eg. block, euclidean)
%         -- sigma
%              Kernel width (k nearest neighbor)
%         -- lambda
%              Regulariztion parameter
%         -- kernel
%              Kernel
%         -- maxit
%              Maximum number of iteration
%         -- folds
%              Number of folds
%         -- soft
%              Kernel type [1: Exponential kernel; 0: limited kernel]
%         -- plot
%              Flag of showing the learning process [1: Plot of the learning
%              process; 0: do not plot]
%--------------------------------------------------------------------------
% OUTPUT:
%     ACC: Crossvalidation accuracy of each fold
%     w_history: Feature weight of each fold
%     Predict_value: Prediction value
%--------------------------------------------------------------------------
% by Lu Li
% update history: 7/23/2018
%% ==========================================================================
% partitioning samples into subsets
fold_id = selectFOLD(label, Para.folds);
nk = length(unique(label)); % number of label class

%% initializations
ACC = zeros(1,max(fold_id));
w_history = zeros(size(X,1),10);
Predict_value = zeros(size(label));
Para.kernel = 'parabolic';

%% apply crossvalidation
for k=1:max(fold_id)
    display(['Forld ' num2str(k)]);
    Para.plotfigure = 0;
    % partitioning samples into training and testing
    train_patterns = X(:,fold_id~=k);
    train_targets = label(fold_id~=k);
    test_patterns = X(:,fold_id==k);
    test_targets = label(fold_id==k);
    Para_train = Para;
    
    %% feature selection using LOGO
    [Weight,~,~] = Logo_kernel(train_patterns, train_targets, Para_train);
    Weight = Weight(:);
    w_history(:,k) = Weight;
    
    %% classification
    index_all = cell(1,nk);
    N = zeros(1,nk);
    patterns_nk = cell(1,nk);
    for i=1:nk
        index_all{i} = find(train_targets == i);
        N(i) = length(index_all{i});
        patterns_nk{i} = train_patterns(:, index_all{i});
    end
    Weight = Weight(:);
    Prho = zeros(size(test_patterns,2),nk);
    
    % predict based on the training data
    for n = 1:size(test_patterns,2)
        test = test_patterns(:,n);
        for i=1:nk
            % calculate the distance to samples with same label
            temp_H = abs(patterns_nk{i}-test*ones(1,N(i)));
            dist_H = (Weight)'*temp_H;
            % calculate the distance to samples with different label
            ID_M = find(train_targets~=i);
            ID_n = length(ID_M);
            patterns_M = train_patterns(:,ID_M);
            temp_M = abs(patterns_M-test*ones(1,ID_n));
            dist_M = (Weight)'*temp_M;
            
            prob_H = kernel_fun(dist_H,Para_train);
            prob_M = kernel_fun(dist_M,Para_train);
            
            rho = sum(dist_M.*prob_M)-sum(dist_H.*prob_H);
            Prho(n,i) = 1/(1+exp(-rho));
        end
    end
    
    % compute testing error
    Prediction = zeros(size(Prho,1),1);
    for i = 1:length(Prediction)
        [~,Prediction(i)] = max(Prho(i,:),[],2);
    end
    Predict_value(fold_id==k) = Prediction;
    Test_Error = length(find(Prediction(:)~=test_targets(:)))/length(test_targets);
    Test_Result = Test_Error;
    R = (Test_Result*100);
    ACC(k) = 100-R;
end
end


