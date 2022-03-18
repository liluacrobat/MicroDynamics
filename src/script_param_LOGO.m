function [opt_lambda, ACC_LOGO] = script_param_LOGO(Table_otu, Table_clinic, params)
%% ========================================================================
% Determine the regularization parameter of LOGO using 10-fold crossvalidation
%
%--------------------------------------------------------------------------
% Input
%   Table_otu    : OTU table after preprocessing
%       -- rel
%            Relative abundance of OTUs
%       -- logRel
%            Relative abundance of OTUs after 10-base log transformation
%       -- tax
%            Taxonomy
%   Table_clinic : Meta table of clinic information
% 
%   params       : Parameters
%       -- lam_ls
%            Range of the regularization parameter lambda [default: 10^-5~100]
%       -- sigma
%            Kernel width (k nearest neighbor) [default: 10]
%       -- folds
%            Number of folds for cross-validation [default: 10]
%--------------------------------------------------------------------------
% Output
%   opt_lambda : optimal value of the
%   ACC_LOGO   : accuracy of cross-validation
%--------------------------------------------------------------------------
% Author: Lu Li
% update history: 08/10/2020
%% ========================================================================
if nargin<3
    params.sigma = 10;
    params.lam_ls =logspace(-5,2,15);
    params.folds = 10;
else
    if ~isfield(params, 'sigma')
        params.sigma = 10;
    end
    if ~isfield(params, 'lam_ls')
        params.lam_ls =logspace(-5,2,15);
    end
    if ~isfield(params, 'folds')
        params.folds = 10;
    end
end

%% initializations
params.sigma = params.sigma+1;
lam_ls = params.lam_ls;
ACC_LOGO = zeros(1,length(lam_ls));
training = Table_otu.logRel;
label = Table_clinic.label;

% remove the OTUs with counts <=20 to facilitate the feature selection
count = Table_otu.count;
idx = sum(count,2)>20;
training = training(idx,:);

%% 10-fold corss-validation
for pi = 1:length(lam_ls)
    para4logo = params;
    para4logo.lambda = lam_ls(pi);
    [ACC{pi}, ~, ~] = script_crossvail_LOGO(training, label, para4logo);
    ACC_LOGO(pi) = mean(ACC{pi});
end
[~, idx_best] = max(ACC_LOGO);
opt_lambda = lam_ls(idx_best);
end