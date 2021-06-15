function [label_sorted, level_sorted] = sort_label(label, label_level, order)
%% ========================================================================
% Sort labels by the given order
%
%% ========================================================================
label = label-min(label)+1;
label_sorted = order(label);
level_sorted = label_level(order);
end