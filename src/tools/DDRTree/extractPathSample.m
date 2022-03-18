function [ordering, onPathDistSample] = extractPathSample(one_path, projectionConnect, Data, Pcurve)
% Extract samples along a path by reordering projected samples according to 
% the ordering of one_path based on principal pcurve 
% 
% Input:
%         one_path: a path based on pcurve, each pt is a curve pt
%projectionConnect: pcurve pt and sample pt connection matrix, each row is
%                   for one pcurve pt, each col is for one sample pt
%                   1: sample pt projected on pcurve pt
%             Data: Sample-by-feature data matrix
%              MST: Minimal spanning tree constructed from pcurve
% output:
%         ordering: ordering of all samples projected onto one_path

num_proj = sum(projectionConnect(one_path, :),2);% No of samples projected 
                                                 % onto one_path of pcurve
% Delete pt w/o projections
one_path(num_proj==0) = []; 
num_proj(num_proj==0) = [];

ordering = zeros(1, sum(num_proj));
pathDist = zeros(1, length(one_path));      % distance of each path pt from the first node
onPathDistSample = zeros(1, sum(num_proj)); % On-path distance of each sample from the first node

% Calculate on-path distance of each path pt
pathDist(1) = 0;
for n = 2:length(one_path)
    pathDist(n) = sqrt(sum((Pcurve(one_path(n),:) - Pcurve(one_path(n-1),:)).^2));% Euclidean distance
    pathDist(n) = pathDist(n) + pathDist(n-1);
end

if(all(num_proj==1))
    % one-2-one correspondence b/w pcurve pt and sample pt    
    for n = 1:length(one_path)
        ordering(n) = find(projectionConnect(one_path(n),:)==1);
    end
    onPathDistSample = pathDist;
else
% multiple sample pt projected onto one pcurve pt
    curr_index = 1;
    for n = 1:length(one_path)
        if(num_proj(n) == 1)
            ordering(curr_index) = find(projectionConnect(one_path(n),:)==1);
            onPathDistSample(curr_index) = pathDist(n);
        else
            projectedSample = find(projectionConnect(one_path(n),:)==1);
            % order the samples projected onto one pcurve pt
            % Find another closest pt one the curve for fitting a line
            if(n==1) % The first point has no predecessor, choose its first successor
                secondIndex = one_path(n+1);
            else
                secondIndex = one_path(n-1);
            end
            [~, theta] = project2line(Pcurve(one_path(n),:)', Pcurve(secondIndex,:)', Data(projectedSample,:)');
            [~, orderIndex] = sort(theta, 'descend');
            ordering(curr_index:curr_index+num_proj(n)-1) = projectedSample(orderIndex);
            onPathDistSample(curr_index:curr_index+num_proj(n)-1) = pathDist(n);
        end
        curr_index = curr_index + num_proj(n);
    end

end


