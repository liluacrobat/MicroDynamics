function [Weight,History,theta] = Logo_kernel(patterns, targets, Para)

%Logo_MulticlassProblem: logo algorithm for feature selction for multiclassication problem
%Y. Sun, S. Todorovic, and S. Goodison,
%Local Learning Based Feature Selection for High Dimensional Data Analysis
%IEEE Trans. on Pattern Analysis and Machine Intelligence, vol. 32, no. 9, pp. 1610-1626, 2010.
%--------------------------------------------------------------------------
%INPUT:
%     patterns:  training data: [x1,x2,...xn] Each column is an observation
%      targets:  class label = {1,2,...,C}
%         Para:  parameters. More specifically,
%   Para.sigma:  kernel width
%  Para.lambda:  regulariztion parameter
%    Para.plotfigure:  1: plot of the learning process; 0: do not plot
%OUTPUT:
%       Weight:  weight of features
%--------------------------------------------------------------------------
%by Lu Li, Yijun Sun @University at Buffalo
%update history: July 23, 2018
%% ==========================================================================
Para.kernel = 'parabolic';
Para.maxit = 10;
Para.distance = 'block';
Para.plotfigure = 0;
Para.soft = 0;
distance = Para.distance; % distance metric
sigma = Para.sigma;                 % number of nearest neighbors
lambda = Para.lambda;               % regulariztion parameter
plotfigure = Para.plotfigure;       % whether the progress and feature weights are plotted

Uc = unique(targets);
if min(Uc)~=1 || max(Uc)~=length(Uc)
    error('Targets should run from 1 to C !')
end
[dim,N_patterns] = size(patterns);
N=zeros(1,length(Uc));
for n=1:length(Uc)
    temp = find(targets==n);
    index{n} =temp;
    N(n) = length(temp);
end

Original_dim = dim;
Original_index = 1:dim;
History = [];
Weight =  1/sqrt(dim)*ones(dim,1); %initial guess
History(:,1) = Weight;
P.lambda = lambda;
Difference =1;t=0;theta =[];
index_all = 1:N_patterns;
maxit=Para.maxit;
while  Difference>0.001 && t<=maxit
    
    t = t+1;
    NM = zeros(dim,N_patterns);
    NH = zeros(dim,N_patterns);
    V = (Weight(:).^2)';
    
    for i = 1:N_patterns
        index_SameClass = setdiff(index{targets(i)},i);
        
        index_DiffClass = setdiff(index_all,index{targets(i)});
        
        switch lower(distance)
            case {'euclidean'}
                Temp_SameClass            = (patterns(:,index_SameClass) - patterns(:,i)*ones(1,length(index_SameClass))).^2;
                Temp_DiffClass            = (patterns(:,index_DiffClass) - patterns(:,i)*ones(1,length(index_DiffClass))).^2;
            case {'block'}
                Temp_SameClass            = abs(patterns(:,index_SameClass) - patterns(:,i)*ones(1,length(index_SameClass)));
                Temp_DiffClass            = abs(patterns(:,index_DiffClass) - patterns(:,i)*ones(1,length(index_DiffClass)));
        end
        
        
        if t==1
            dist_SameClass    = sum(Temp_SameClass,1)/sqrt(dim);
            dist_DiffClass    = sum(Temp_DiffClass,1)/sqrt(dim);
        else
            dist_SameClass    = (V)*Temp_SameClass;
            dist_DiffClass    = (V)*Temp_DiffClass;
        end
        
        if Para.soft==1
            % temp_index_SameClass = find(dist_SameClass==0);
            prob_SameClass = exp(-dist_SameClass/sigma);%prob_SameClass(temp_index_SameClass) = 0;
            if sum(prob_SameClass)~=0;prob_S = prob_SameClass/sum(prob_SameClass);else;[dum,I] = sort(dist_SameClass);prob_S=zeros(size(I));prob_S(I(1))=1;end
            prob_DiffClass = exp(-dist_DiffClass/sigma);
            if sum(prob_DiffClass)~=0;prob_D = prob_DiffClass/sum(prob_DiffClass);else;[dum,I] = sort(dist_DiffClass);prob_D=zeros(size(I));prob_D(I(1))=1;end
        else
            % Assign the probability of NM/NH with k nearest neighbors
            prob_S=kernel_fun(dist_SameClass,Para);
            prob_D=kernel_fun(dist_DiffClass,Para);
        end
        NH(:,i) = Temp_SameClass*prob_S(:);
        NM(:,i) = Temp_DiffClass*prob_D(:);
    end
    
    Z = NM-NH;
    CostDiff = 1000; Cost(1) = 10000;
    j=1;
    while CostDiff>0.01*Cost(j) && length(Cost)<500
        j= j+1;
        a = (Weight.^2)'*Z; % Margin
        Result = 1./(1+exp(a));
        descent = lambda*Weight-(Z*Result(:)).*Weight;
        
        P.Weight = Weight;
        P.descent = descent;
        [alpha, Cost(j)] = fminbnd(@(p) logocost_multiclass_cost(p, P, Z), 0, 1);
        
        Weight = Weight-alpha*descent;
        CostDiff = abs(Cost(j)-Cost(j-1));
    end
    
    Weight = abs(Weight);
    Difference = norm(abs(Weight/max(Weight)-History(:,t)/max(History(:,t))));%max(abs(Weight/max(Weight)-History(:,t)/max(History(:,t))));
    theta(t) = Difference;
    beta(t)=Cost(j);
    History(:,t+1) = Weight;
    
    if t==1;index_zeros = [];end%find(Weight<=10^(-5));end
    if t>=2;index_zeros = [];end%find(Weight<=10^(-5));end
    patterns(index_zeros,:)=[];
    dim = size(patterns,1);
    Weight(index_zeros)=[];
    History(index_zeros,:)=[];
    Original_index(index_zeros)=[];
end
temp = zeros(1,Original_dim);
temp(Original_index) = Weight.^2;
Weight = temp;

%Monitoring the feature weights
if plotfigure ==1
    figure;
    semilogy(theta,'-o','LineWidth',1,'MarkerFaceColor','w','MarkerSize',10)
    title('Theta');
    xlabel('Number of Iterations');
    ylabel('Difference')
    grid on
    drawnow
end

return
%% ==================End of the code===================================
end
