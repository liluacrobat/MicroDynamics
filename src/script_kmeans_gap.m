function [cluster_num, eva_Gap] = script_kmeans_gap(Feature_Table)
%% ========================================================================
% Estimate the number of clusters based on gap statistics using k-means
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
%--------------------------------------------------------------------------
% Output
%   cluster_num  : Optimal number of clusters
%   eva_Gap      : Gap statistics
%--------------------------------------------------------------------------
% Author: Lu Li
% update history: 08/10/2020
%% ========================================================================
%% clustering analysis using K-means
% determine the number of clusters using gap statistics

X = Feature_Table.logRel;
[cluster_num, eva_Gap] = gapKmeans(X);

% Plot Gap results
plotGap(cluster_num, eva_Gap);
pbaspect([2 1 1]);
end

function [h, x, y] = plotGap(cluster_num, eva_Gap)
% Plot Gap statistic result

Gap = eva_Gap.ExpectedLogW-eva_Gap.LogW;
Gap_SE = Gap-eva_Gap.SE;
delta_Gap = Gap(1:end-1)-Gap_SE(2:end);

% Data to plot
x = eva_Gap.InspectedK(1:end-1);
y = delta_Gap;
x = x(2:end);
y = y(2:end);
h = figure; hold on;
bar(x, y);
xlabel('Number of Clusters');
ylabel('Gap(k)-Gap(k+1)+SE(k+1)');
set(gca,'FontSize',14)
pbaspect([2.8 1 1])
end

function [numClusters, eva_Gap] = gapKmeans(Data)
%============================================================%
% Perform kmeans clustering on data after feature selection
%============================================================%

%% ======================================%
% List of options and parameters
%======================================%
% @@ Options @@
para.clusterAlg = 'Kmeans'; % use consensus clustering
% @@ Parameters @@
para.fs_threshold = 1e-2;
para.CLUSTER_NUM_CHOICES = 1:10; % Candidate numbers of clusters
%% =====================================%
% Load and prepare data
%=======================================%
DATA = Data;
% load annotation

%% ======================================%
% Clustering
% =======================================%
%Determine the number of clusters using Gap statistics
% parpool;
rng(25);
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton', 'replicate',20, ...
    'Options',statset('UseParallel',1)));

tic;
eva_Gap = evalclusters(transpose(DATA),'kmeans','gap','KList',para.CLUSTER_NUM_CHOICES, ...
    'SearchMethod', 'firstMaxSE', 'ReferenceDistribution','PCA', 'B',100); %%PCA

toc;
%%
numClusters = eva_Gap.OptimalK;
display(['Optimal # of Clusters: ' num2str(numClusters)]);
end




