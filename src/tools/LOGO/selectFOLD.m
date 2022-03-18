function fold_id = selectFOLD(label, folds)
%% ==========================================================================
%
% selectFOLD: partitioning samples into subsets
%
%--------------------------------------------------------------------------
%INPUT:
%     label:  class label = {1,2,...,C}
%     seed: random seed
%     folds: number of folds fro crossvalidation
%
%OUTPUT:
%     fold_id: partition of folds
%
%--------------------------------------------------------------------------
% by Lu Li
% update history: 7/23/2018
%% ==========================================================================
idx = 1:length(label);
rng(98,'twister');
nQ = length(unique(label));
l_G = zeros(1,nQ);
perm_G = cell(1,nQ);
step_G = zeros(1,nQ);
idx_G = cell(1,nQ);
% random select samples for each class
for i = 1:nQ
    l_G(i) = length(find(label==i));
    perm_G{i} = randperm(l_G(i));
    step_G(i) = round(l_G(i)/folds);
    idx_G{i} = idx(label==i);
end
% assign samples for each fold
fold_id = ones(length(label),1)*folds;
for i = 1:(folds-1)
    for j = 1:length(unique(label))
        idx_B = idx_G{j};
        perm_B = perm_G{j};
        l_B = l_G(j);
        fold_id(idx_B(perm_B(round(l_B/folds*(i-1))+1:round(l_B/folds*i)))) = i;
    end
end
%% ==================End of the code===================================
end