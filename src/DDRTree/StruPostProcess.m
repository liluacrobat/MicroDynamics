function [extracted_path, extracted_pathDist,Pcurve_etd,extracted_curve,...
    edges_etd,DataProjection] = StruPostProcess(Pcurve,Data,subtype)
% Main function of post-processing the principal tree
%
% Input:
%      Pcurve: N-by-Feature matrix where N is # of points on the principal
%      tree
%      Data: M-by-Feature matrix where M is # of data samples where
%            pcurve is learned
%      subtype: subtype label of samples
% Output:
%      extracted_path: extracted path including the ordered index of
%      samples
%      extracted_pathDist: extracted path including the progression
%      extracted_curve: extracted path including the ordered index of
%      principal tree points
%      distance of each samples
%      DataProjection: projection of data on the principalk tree
%      Pcurve_etd:  extended principal curve
%% remove extremely small branches
normal_label = 'HC';
Pcurve = Pcurve';
Data = Data';
Pcurve = removeSmallBranch(Pcurve, Data, subtype, normal_label);
if size(Pcurve,1)>1
    %% ====================================================
    % Extend the curves
    % =====================================================
    MST = buildMST(Pcurve);
    % search for the root within HC
    root_index = search4Root(Pcurve, Data, subtype, normal_label);
    
    Pcurve = extendPcurve(Pcurve, Data, root_index(1), subtype, normal_label);
    % search for the root after extending curve
    root_index = search4Root(Pcurve, Data, subtype, normal_label);
    
    %build a MST from pcurve
    MST = buildMST(Pcurve);
    degree = sum(full(MST)>0, 1);
    display(['Number of end points:' num2str(length(find(degree==1)))]);
    reassign = true;
    if(reassign)
        % reassign samples to its nearest neighbors in the curve
        distanceMatrix = pdist2(Data, Pcurve);
        projectionConnect = zeros(size(Pcurve, 1), ...% row: pcurve pt, col: data pt
            size(Data, 1));     % 1: data pt projected on pcurve pt
        for n = 1:size(Data, 1)
            [~, index] = min(distanceMatrix(n,:));
            DataProjection(n,:) = Pcurve(index,:);
            projectionConnect(index,n) = 1;
        end
    end
    
    %%
    % build a MST from pcurve
    MST = buildMST(Pcurve);
    degree = sum(full(MST)>0, 1);
    singleton = find(degree==1);
    singleton = setdiff(singleton,root_index);
    num_path = length(singleton);
    num_path_all = num_path*length(root_index);
    %% ====================================================
    % Extract progression paths from root to every end-point
    % ====================================================
    extracted_path = cell(length(root_index), num_path);
    extracted_pathDist = cell(length(root_index), num_path);
    extracted_curve = cell(length(root_index), num_path);
    num_path = length(singleton);
    for k=1:length(root_index)
        for n = 1:num_path
            [~,  one_path, ~] = graphshortestpath(MST, root_index(k), singleton(n));
            extracted_curve{k,n} = one_path;
            [extracted_path{k,n}, extracted_pathDist{k,n}] = extractPathSample(one_path, projectionConnect, Data, Pcurve);
        end
        path_names(k,:) = cellfun(@num2str, num2cell((1:num_path)+(k-1)*num_path), 'UniformOutput', false);
    end
    Pcurve_etd = Pcurve';
    edges_etd = MST;
else
    disp('Error');
end
end
function exPcurve= extendPcurve(Pcurve, Data, root_index, subtype, normal_label)
% Extend pcurve to resovle the problem of many sample pts
% mapped to the same curve end point
% Input:
%      Pcurve: N-by-Feature matrix where N is # of curve pts
%      Data: M-by-Feature matrix where M is # of data samples where
%            pcurve is learned
%      root_index: Root node index on Pcurve as a path staring pt
%      subtype: subtype label of each sample
%      normal_label: the label of normal
% Output:
%      exPcurve: Extended curve pts
%DataProjection: projected data pt on the extended curve
%projectionDist: Distance b/w sample pt with its projection pt
%eXprojectionConnect: projectionConnect matrix after curve extension

small_branch = 0.1; % branch_path_ratio smaller than less value is a small branch
fit_sample   = 0.7;   % Percentage of sample on a branch used to extend curve


DataProjection = zeros(size(Data));
exPcurve = Pcurve;
distanceMatrix = pdist2(Data, Pcurve);
projectionDistance = zeros(1, size(Data, 1));
projectionConnect = zeros(size(Pcurve, 1), ...% row: pcurve pt, col: data pt
    size(Data, 1));     % 1: data pt projected on pcurve pt

% Find the nearest curve pt as sample's projection pt
for n = 1:size(Data, 1)
    [dist, index] = min(distanceMatrix(n,:));
    projectionDistance(n) = dist;
    DataProjection(n,:) = Pcurve(index,:);
    projectionConnect(index,n) = 1;
end
eXprojectionConnect = projectionConnect; % projectionConnect after exten

% Build a MST from pcurve
MST = buildMST(Pcurve);

% Find points to extend
extend_end_point = true;
if(extend_end_point)
    degree = sum(full(MST)>0, 1);
    singleton = find(degree==1);
else
    num_proj = sum(projectionConnect,2);
    [~, index] = sort(num_proj, 'descend');
    index(num_proj(index) <= 20) = []; % Remove end points with small amount of pts
    singleton = index;
end
singleton=setdiff(singleton,root_index);
% Plot curve and end points to be extended from
plot_curve = 0;
if plot_curve == 1
    figure; hold on; grid on;
    plot3(Pcurve(:,1), Pcurve(:,2), Pcurve(:,3), 'o', ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w');
    plot3(Pcurve(singleton,1), Pcurve(singleton,2), Pcurve(singleton,3), 'pr', 'MarkerSize', 20);
    set(gca, 'ZColor',[1 1 0],'YColor',[1 1 0],'XColor',[1 1 0],'Color',[0 0 0]);
    for n = 1:length(singleton)
        sampleIndex = find(projectionConnect(singleton(n),:)==1);
        plot3(Data(sampleIndex,1), Data(sampleIndex,2), Data(sampleIndex,3), 'oc');
    end
    title('Extend the principal curve','FontSize',18)
    % Find all branches in MST
end
[branch, branch_ratio] = search4Branch(MST, root_index, singleton);
if length(branch)==1
    fit_sample=0.3;
end
[branchMain, branch_ratioMain] = search4MainPath(MST, root_index, singleton);
singleton=[singleton root_index];
branch=[branch branchMain];
branch_ratio=[branch_ratio branch_ratioMain];
% Extend curve

for n = 1:length(singleton)
    
    sampleIndex = find(projectionConnect(singleton(n),:)==1);
    
    sub_branch = Pcurve(branch{n},:);
    if branch_ratio(n)>0.05
        if plot_curve == 1
            plot3(sub_branch(:,1), sub_branch(:,2), sub_branch(:,3), '-ob', 'LineWidth',2, 'MarkerSize', 4);
        end
        poly_degree = 3;
        if(branch_ratio(n) > small_branch)
            if(fit_sample < 1)
                sub_branch = sub_branch(round(length(branch{n})*(1-fit_sample)):end,:);
            end
        else
            % Small branch fit a line
            poly_degree = 1;
        end
        if plot_curve == 1
            plot3(sub_branch(:,1), sub_branch(:,2), sub_branch(:,3), '-oy', 'LineWidth',2, 'MarkerSize',6);
        end
        if ~isempty(sampleIndex)
            [projection, extended, extendConnect] = project2poly(sub_branch, Data(sampleIndex,:), poly_degree);
            
            % Extended curve
            exPcurve = [exPcurve; extended];
            if plot_curve == 1
                plot3(extended(:,1), extended(:,2), extended(:,3), 'yo', 'MarkerSize', 6);
            end
            % Assign projection
            DataProjection(sampleIndex,:) = projection;
            if plot_curve == 1
                plot3(projection(:,1), projection(:,2), projection(:,3), ...
                    'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'MarkerSize', 6);
            end
            %recalculate  projection distance
            projectionDistance(sampleIndex) = sqrt(sum((Data(sampleIndex,:)-projection).^2, 2));
            
            % Re-assign projectionConnection
            eXprojectionConnect(singleton(n),sampleIndex) = 0;
            tempConnect = zeros(size(extended, 1), size(Data, 1));
            tempConnect(:, sampleIndex) = extendConnect;
            eXprojectionConnect = [eXprojectionConnect; tempConnect]; % Append to end
        end
    end
end
% Remove the small branches generated by extending the curve
exPcurve = removeSmallBranch(exPcurve, Data, subtype, normal_label);
end

function [newPcurve,root]= removeSmallBranch(Pcurve, Data, subtype, normal_label)
%% ====================================================
%  Remove extremely small branches
%% ====================================================
distanceMatrix = pdist2(Data, Pcurve);
projectionConnect = zeros(size(Pcurve, 1), ...% row: pcurve pt, col: data pt
    size(Data, 1));     % 1: data pt projected on pcurve pt

% Find the nearest curve pt as sample's projection pt
for n = 1:size(Data, 1)
    [~, index] = min(distanceMatrix(n,:));
    projectionConnect(index,n) = 1;
end

% Search for root node
root_index = search4Root(Pcurve, Data, subtype, normal_label);
% Build a MST from pcurve
MST = buildMST(Pcurve);
branch2remove = [];
% Find all branches in MST
[all_branch, branch_ratio, branch_length] = search4Branch(MST, root_index(1));
if length(all_branch)>1
    % Rmove branches with small branches
    for n = 1:length(all_branch)
        branch = all_branch{n};
        if branch_ratio(n)< 0.03%|branch_length(n)<0.5
            branch2remove = [branch2remove, branch(2:end)];
        end
    end
    newPcurve = Pcurve;
    newPcurve(branch2remove,:) = [];
else
    newPcurve = Pcurve;
end
end
function singleton=FindSingleton(Pcurve)
MST = buildMST(Pcurve);
degree = sum(full(MST)>0, 1);
singleton = find(degree==1);
end
function plotTree(X,F,Y,stree,root)
plotFIGURE([X F],[Y;ones(size(F,2),1)*(max(Y)+1)]);
U=sort(unique(Y));
str=cell(1,length(U));
FaceColor =  distinguishable_colors(length(U));
mapped_data=X;
figure,
for i=1:length(U)
    plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:));
    hold on
    str{i}=num2str(U(i));
end

% hold on
%     plot(F(1,:),F(2,:),'or','MarkerFaceColor','r');
hold on
[m,n]=size(stree);
for i=1:m
    for j=1:n
        if stree(i,j)~=0
            plot3([F(1,i) F(1,j)],[F(2,i) F(2,j)],[F(3,i) F(3,j)],'-r','linewidth',4);
        end
    end
end
plot3(F(1,root),F(2,root),F(3,root),'pr','MarkerSize',16);
hold off
grid
end
function plotFIGURE(X,Y)
if nargin<2
    Y=ones(size(X,2),1);
end
U=sort(unique(Y));
str=cell(1,length(U));
FaceColor =  distinguishable_colors(length(U));
D=size(X,1);
mapped_data=X;
switch D
    case 2
        figure,
        for i=1:length(U)
            plot(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:),'linewidth',3);
            hold on
            str{i}=num2str(U(i));
        end
        legend(str)
        grid
    case 3
        figure,
        for i=1:length(U)
            plot3(mapped_data(1,Y==U(i)),mapped_data(2,Y==U(i)),mapped_data(3,Y==U(i)),'o','MarkerFaceColor',FaceColor(i,:));
            hold on
            str{i}=num2str(U(i));
        end
        legend(str)
        grid
end
end
function [branch, branch_ratio,branch_length] = search4Branch(MST, root_index, singleton)
% Find all branches from a minimal spanning tree (MST)
% Input:
%       MST: a sparse connection matrix to represent a MST
% Output:
%    branch: pt index of all branches found in MST

degree = sum(full(MST)>0, 1);
% No specified target points, use end points
if(nargin < 3)
    singleton = find(degree==1);
    singleton = setdiff(singleton,root_index);
end

branch = cell(1, length(singleton));
branch_length = zeros(1, length(singleton));
path_length = zeros(1, length(singleton));
% Find branch each end point resides
for n = 1:length(singleton) % Traverse every path from root node
    [path_length(n),  one_path, ~] = graphshortestpath(MST, root_index, singleton(n));
    
    
    % calculate branch
    branch_degree = degree(one_path);
    cross_index = find(branch_degree>2); % cross plot
    if(isempty(cross_index))
        % No cross point on this path, set branch length as one path length
        branch{n} = one_path;
        branch_length(n) = path_length(n);
        %         if(n==2)
        %            branch{1} = fliplr(one_path); % for root node
        %         end
    else
        % No. of samples from end point to the nearest cross point
        branch{n} = one_path(cross_index(end):end);
        [branch_length(n),  ~, ~] = graphshortestpath(MST, one_path(cross_index(end)), one_path(end));
        %           if(n==2)
        %               branch{1} = fliplr(one_path(1:cross_index(1))); % for root node
        %           end
    end
end
branch_ratio = branch_length./path_length;

end
function [branch, branch_ratio] = search4MainPath(MST, root_index, singleton)
% Find all branches from a minimal spanning tree (MST)
% Input:
%       MST: a sparse connection matrix to represent a MST
% Output:
%    branch: pt index of all branches found in MST

degree = sum(full(MST)>0, 1);
% No specified target points, use end points
if(nargin < 3)
    singleton = find(degree==1);
end

branch = cell(1, length(singleton));
branch_length = zeros(1, length(singleton));
path_length = zeros(1, length(singleton));
% Find branch each end point resides
flag=0;
for n = 1:length(singleton) % Traverse every path from root node
    [path_length(n),  one_path, ~] = graphshortestpath(MST, root_index, singleton(n));
    
    
    % calculate branch
    branch_degree = degree(one_path);
    cross_index = find(branch_degree>2); % cross plot
    if(~isempty(cross_index))
        new_root=one_path(cross_index(1));
        [branch, branch_ratio] = search4Branch(MST, new_root, root_index);
        flag=1;
        break;
    end
end
if flag==0
    new_root=one_path(end);
    [branch, branch_ratio] = search4Branch(MST, new_root, root_index);
end
end
