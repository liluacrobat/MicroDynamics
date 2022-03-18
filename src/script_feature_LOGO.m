function Feature_Table = script_feature_LOGO(Table_otu, Table_clinic, params)
%% ========================================================================
% Feature selection within LOGO framework
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
%       -- sigma
%            Kernel width (k nearest neighbor) [default: 10]
%       -- lambda
%            Regularization parameter [default: 10^-4]
%       -- threshold
%            Threshold of feature weight
%--------------------------------------------------------------------------
% Output
%   Feature_Table : Table of the selected OTUs
%       -- tax  
%            Taxonomy of the selected OTUs
%       -- logRel 
%            Relative abundance of the selected OTUs after 10-base log 
%            transformation
%       -- weight
%            Feature weights of the selected OTUs learned from LOGO
%       -- weight_raw
%            Feature weights of all the OTUs learned from LOGO
%--------------------------------------------------------------------------
% Author: Lu Li
% update history: 08/10/2020
%% ========================================================================

if nargin<3
    params.sigma = 10;
    params.lam_ls =logspace(-5,2,15);
    params.threshold = 0.001;
else
    if ~isfield(params, 'sigma')
        params.sigma = 10;
    end
    if ~isfield(params, 'lam_ls')
        params.lam_ls = 10^-4;
    end
    if ~isfield(params, 'threshold')
        params.threshold = 0.001;
    end
end

%% initializations
params.plotfigure = 0;
params.sigma = params.sigma+1;
training = Table_otu.logRel;
label = Table_clinic.label;

% remove the OTUs with counts <=20 to facilitate the feature selection
count = Table_otu.count;
Weight_raw = zeros(size(count,1),1);
idx = sum(count,2)>20;
training = training(idx,:);
tax = Table_otu.tax(idx,:);

% feature selection using LOGO
[Weight,~,~] = Logo_kernel(training, label, params);
Weight_raw(idx) =  Weight;
Feature_Table.weight_raw = Weight_raw;

% normalize the feature weight
Weight = Weight/sum(Weight);

ids = Weight>params.threshold;
Feature_Table.tax = tax(ids);
Feature_Table.logRel = training(ids,:);
Feature_Table.weight = Weight(ids);

[Weight_sorted,~]=sort(Weight,'descend');
figure,semilogx(Weight_sorted,'-k','linewidth',2);
hold on
plot([sum(ids) sum(ids)],[0 max(Weight)],'--r','linewidth',2);
plot([1 sum(ids)],[0.001 0.001],'--r','linewidth',2);
ylabel('Normalized Feature Weight');
xlabel('OTU Index');
set(gca,'FontSize',14);
axis([1 10^5 0 max(Weight)])
end


