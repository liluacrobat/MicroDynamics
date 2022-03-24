function cidx = script_consensus_clustering(Feature_Table, params)
%% ========================================================================
% Random sampling based consensus clustering
%
%--------------------------------------------------------------------------
% Input
%   Feature_Table : Table of the selected OTUs
%       -- tax
%            Taxonomy of the selected OTUs
%       -- logRel
%            Relative abundance of the selected OTUs after 10-base log
%            transformation
%       -- weight
%            Feature weight of the selected OTUs learned from LOGO
%   params        : Parameters
%       -- med          : Method used for clustering 
%                         options 
%                            km: k-means
%                            dmm: dirichlet multinomial mixtures
%                            cst: community state typing with partition 
%                                 around medoids
%       -- alpha        : Parameter for Dirichlet process prior
%       -- cluster_num  : Number of clusters (available when med=km)
%       -- iters        : Number of iterations for consensus clustering
%--------------------------------------------------------------------------
% Output
%   cidx : Cluster labels
%--------------------------------------------------------------------------
% Author: Lu Li
% update history: 03/18/2022
%% ========================================================================

cluster_num = params.cluster_num;
X = Feature_Table.logRel;
[consensus_mtx,opt] = ConsensusClustering(X, params);
switch params.med
    case 'dmm'
        cluster_num = opt;
end
% perform hierarchical clustering on the consensus matrix
Y = pdist(1-consensus_mtx);
Z = linkage(Y,'complete');
cidx = cluster(Z,'maxclust', cluster_num);
end

function [consensus,opt] = ConsensusClustering(X,params)
opt = [];
if isfield(params,'iters')
    iter = params.iters;
else
    iter = 1000;
end

switch params.med
    case 'km'
        K = params.cluster_num;
    case 'dmm'
        nc = zeros(1,iter);
    case 'cst'
        K = params.cluster_num;
end

X=X';
rng(11,'twister');
[m,~] = size(X);
ms = round(m*0.8);
CMmat = zeros(m);
CMCmat = zeros(m);

for i=1:iter
    idx = randperm(m,ms);
    Xi = X(idx,:);
    CMmat(idx,idx) = CMmat(idx,idx)+1;
    switch params.med
        case 'km'
            [cidxi, ~] = kmeans(Xi, K, 'Replicates',20,'OnlinePhase','off');
        case 'dmm'
            para4dp.alpha = params.alpha;
            [cidxi, ~, ~, ~] = mixGaussGb(Xi',para4dp);
            nc(i) = max(cidxi);
        case 'cst'
            cidxi = kmedoids(Xi, K);
    end
    
    temp = zeros(ms);
    for j=1:ms
        sel = cidxi==cidxi(j);
        temp(j,sel) = temp(j,sel)+1;
    end
    CMCmat(idx,idx) = CMCmat(idx,idx)+temp;
end
switch params.med
    case 'dmm'
        uc = unique(nc);
        nc_sta = zeros(1,length(uc));
        for i=1:length(uc)
            nc_sta(i) = sum(nc==uc(i));
        end
        [~,ids] = max(nc_sta);
        opt = uc(ids);  
end
consensus = CMCmat./CMmat;
end




