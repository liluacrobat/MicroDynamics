function [Table_otu, Table_clinic] = script_data_processing(filen_otu, file_meta, params)
%% ========================================================================
% Load the OTU table and meta data. Exclude samples without enough sequencing
% depth. Filter out OTUs with low counts.
%
%--------------------------------------------------------------------------
% Input
%   filename_otu  : Path to the file of a DxN OTU table, where D is the number of
%                   OTUs and N is the number of samples
%   filename_meta : Path to the meta file of clinic information, where the first
%                   column is the sample ID
% 
%   params        : Parameters
%       -- min_count
%            Number of observation (sequence) count to apply as the minimum
%            total observation count of a sample for that sample to be retained.
%            If you want to include samples with sequencing depth higher than
%            or equal to 10,000, you specify 10,000 [default: 10,000]
%       -- pseudo_count
%            A small number added to the relative abundance before 10-base log
%            transformation [default: 10^-6]
%       -- last_tax
%            Flag of whether the last column of the OTU table is taxonomy
%            or not. If the last column of the table is the taxonomy, you
%            specify 1 [default: 0]
%       -- col_label
%            The Column of the clinical information used for feature
%            selection
%       -- mapping
%            Mapping from class categories to numerical labels
%--------------------------------------------------------------------------
% Output
%   Table_otu : OTU table
%       -- rel
%            Relative abundance of OTUs
%       -- logRel
%            Relative abundance of OTUs after 10-base log transformation
%       -- tax
%            Taxonomy
%   Table_clinic : Meta table of clinic information
%--------------------------------------------------------------------------
% Author: Lu Li
% update history: 08/10/2020
%% ========================================================================

%% Initialization
if nargin<3
    params.min_count = 10000;
    params.pseudo_count = 10^-6;
    params.last_tax = 0;
else
    if ~isfield(params, 'min_count')
        params.min_count = 10000;
    end
    if ~isfield(params, 'pseudo_count')
        params.pseudo_count = 10^-6;
    end
    if ~isfield(params, 'last_tax')
        params.last_tax = 0;
    end
end

%% Load the OTU table and meta data
Table_otu_raw = readtable(filen_otu, 'Delimiter', '\t');
Table_clinic_raw = readtable(file_meta, 'Delimiter', '\t');

var_name = Table_otu_raw.Properties.VariableNames;
Table_otu.sample = var_name(2:end-params.last_tax);

Table_otu.count = table2array(Table_otu_raw(:,2:end-params.last_tax));

if params.last_tax==1
    Table_otu.tax = table2cell(Table_otu_raw(:,end));
else
    Table_otu.tax = table2cell(Table_otu_raw(:,1));
end

%% Filtering based on observation count

counts_per_sample = sum(Table_otu.count,1);
idx_sample = counts_per_sample >= params.min_count;

Table_otu.sample = Table_otu.sample(idx_sample);
Table_otu.count = Table_otu.count(:,idx_sample);
counts_per_sample = counts_per_sample(idx_sample);

idx_OTU = sum(Table_otu.count,2) > 0;
Table_otu.tax = Table_otu.tax(idx_OTU);
Table_otu.count = Table_otu.count(idx_OTU,:);

% Calculate relative abundance and apply 10-base log transformation
[D, N] = size(Table_otu.count);
Table_otu.rel = Table_otu.count./repmat(counts_per_sample, D, 1);

Table_otu.logRel = log10(Table_otu.rel+10^(-6));

Table_clinic = Table_clinic_raw(idx_sample,:);
label_category = table2cell(Table_clinic(:, params.col_label));
n_category = size(params.mapping,1);
label = zeros(N,1);
for i=1:N
    for j=1:n_category
        if strcmpi(label_category{i},params.mapping{j,1})==1
            label(i,1) = params.mapping{j,2};
            break;
        end
    end
end
Table_clinic.label = label;
end