function Demo_MicroDynamics
%% ========================================================================
% A demo of exploring the microbial dynamics in a data set of Crohn's disease
% Lu Li
% update history: 08/10/2020
%% ========================================================================

clc; clear; close all;

%% Setup
addpath(genpath('../src/'));
dataFolder = './data/';

%% 1. Preprocessing
% Load tables of the data set
% Exclude samples without enough sequencing depth
% Exclude OTUs with low counts

file_otu = [dataFolder 'CD_16S_OTU.txt'];
file_meta = [dataFolder 'CD_meta_clinic.txt'];
para4pre.min_count = 10000;
para4pre.pseudo_count = 10^-6;
para4pre.last_tax = 0;
para4pre.col_label = 'Disease_Behavior';

% Mapping from class categories to numerical labels
mapping = {'HC', 1
    'Inflammatory (B1)', 2
    'Stricturing (B2)', 3
    'Penetrating (B3)', 3};
para4pre.mapping = mapping;
[Table_otu, Table_clinic] = script_data_processing(file_otu, file_meta, para4pre);

%% 2. Feature selection
% Use disease behavior to select OTUs associated with CD within LOGO
% framework

para4LOGO.sigma = 10;
para4LOGO.lam_ls = logspace(-5,2,15);
para4LOGO.folds = 10;

% % Optimize the regularization parameter using 10-fold cross-validation
% [opt_lambda, ACC_LOGO] = script_param_LOGO(Table_otu, Table_clinic, para4LOGO);

% load precalculated optimal value of regularization parameter lambda 
load('precalculated/precalculated_step2_feature_selection_cross_validation',...
    'opt_lambda');

para4LOGO.lambda = opt_lambda;
Feature_Table = script_feature_LOGO(Table_otu, Table_clinic, para4LOGO);

%% 3. Random sampling based consensus clustering
% % Estimate the number of clusters based on gap statistics
% [cluster_num, eva_Gap] = script_kmeans_gap(Feature_Table);
% axis([1.1 9.9 -0.024 0.02]);

para4cluster.iters = 1000;
para4cluster.cluster_num = 5;
  
cidx = script_consensus_clustering(Feature_Table, para4cluster);

% reorder the cluster labels
cidx = reorder_cluster(cidx);
save('precalculated/precalculated_step3_clustering','cidx');

%% 4. Embedded structuring learning
% % Optimize the parameters employed in DDRTree using elbow method
% % fix lambda=1, estimate optimal sigma=0.5
% para4ddr.f_sig = 1;
% para4ddr.lambda = 1;
% [sigma_opt, CurveLength_sig, Error_mse_sig] = script_elbow_DDRTree(Feature_Table, para4ddr);
% 
% % fix sigma=0.5, estimate optimal lambda=150
% para4ddr.f_sig = 0;
% para4ddr.sigma = 0.5;
% [lambda_opt, CurveLength_lambda, Error_mse_lambda] = script_elbow_DDRTree(Feature_Table, para4ddr);

para4ddr.sigma = 0.5;
para4ddr.lambda = 150;
para4ddr.col_label = 'Disease_Behavior';
PrincipalTree = script_structure_learning(Feature_Table, Table_clinic, para4ddr);

%% 5. Visualization of results
% plot the principal tree colorcoded by disease behavior
para4visual.FaceColor = [0.2000    0.6275    0.1725
    0.8588    0.4275         0
    0.7137    0.8588    1.0000
    1.0000    0.7137    0.8588];
para4visual.order = [3, 2, 1, 4];
para4visual.post_flag = 1;
para4visual.prog_flag = 1;

annotations = table2cell(Table_clinic(:, 'Disease_Behavior'));
[label, levels]= grp2idx(annotations); % convert annotations into numeric labels

script_visualization(PrincipalTree, label, levels, para4visual);
view(-122,14);

%% ==================End of the code=======================================
end

function cidx = reorder_cluster(cidx)
temp =cidx;
cidx(temp==3) = 1;
cidx(temp==2) = 2;
cidx(temp==5) = 3;
cidx(temp==1) = 4;
cidx(temp==4) = 5;
end