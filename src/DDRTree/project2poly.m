function [projection, branch_extend, projectionConnect, branch_fitted] = project2poly(branch, Data, poly_degree, line_ratio)
% Project Data onto a polynomial line fitted by pcurve in every dimension
% Input:
%       branch: N-by-D curve branch pt
%         Data: M-by-D data sample pt
%  poly_degree: degree of polynomial
%   line_ratio: % qunatile of lambda to estimate extension ratio
% Output:
%   projection: M-by-D data projection on polynomial
%branch_extend: M-by-D extended branch pt
% branch_fitte: M-by-D fitted branch pt
%projectionConnect: M-by-M projection connection matrix

if(nargin < 4)
    line_ratio = 1;
end

% compute branch length
L = transpose(pdist2(branch(1,:), branch));

% Estimate extension ratio
[~, lambda] = project2line(branch(1,:)', branch(end,:)', Data');
extend_ratio = quantile(lambda, line_ratio); % robust to outliers


% Curve fitting of each dimension against length
XX = linspace(L(end)+0.001, max(L)*extend_ratio, size(Data,1));% start from next point
ZZ = linspace(0, max(L)*extend_ratio, size(Data,1)*2);% whole length
branch_extend = zeros(size(Data,1), size(Data,2));
branch_fitted = zeros(size(Data,1)*2, size(Data,2));
% XX = linspace(0, max(L)*extend_ratio, size(Data,1)+size(branch,1));% Extend curve length
% branch_extend = zeros(size(Data,1)+size(branch,1), size(Data,2));

for n = 1:size(branch,2)
    P = polyfit(L, branch(:,n), poly_degree);%3th polynomial
    YY = polyval(P, XX);
    branch_extend(:,n) = YY;
    branch_fitted(:,n) = polyval(P, ZZ);
end


% Calculate projection pt from extended branch
projection = zeros(size(Data));
distanceMatrix = pdist2(branch_extend, Data);
projectionConnect = zeros(size(branch_extend, 1), ...% row:extend pcurve pt col: data
                          size(Data, 1));     % 1: data pt projected on pcurve pt
for n = 1:size(Data, 1)
    [~, index] = min(distanceMatrix(:,n));
    projection(n,:) = branch_extend(index,:);
    projectionConnect(index,n) = 1;
end

% Shrink extended curve to the further pt with data projection
with_proj = find(sum(projectionConnect, 2)>0);
last_index = with_proj(end);
branch_extend(last_index+1:end,:) = [];
projectionConnect(last_index+1:end,:) = [];





 