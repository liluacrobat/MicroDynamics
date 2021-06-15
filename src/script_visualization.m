function script_visualization(PrincipalTree, annotations, levels, params)
%% ========================================================================
% Visualize the principal tree learned from data
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
%   annotations   : Annotation of sampls
%   levels        : Annotation levels
%
%   params        : Parameters
%       -- FaceColor
%            Facecolor for annotations
%       -- order
%            Order used to sort labels
%       -- post_flag
%            Flag of post processing [default: 0]
%              0 : principal tree before post-processing
%              1 : principal tree after post-processing
%       -- prog_flag
%            Plot samples of the extracted progression paths [default: 0]
%--------------------------------------------------------------------------
% Author: Lu Li
% update history: 08/10/2020
%% ========================================================================

%% Initializations
if ~isfield(params, 'post_flag')
    params.post_flag = 0;
end
if ~isfield(params, 'prog_flag')
    params.prog_flag = 0;
end
if params.post_flag==0
    projection = PrincipalTree.projection;
    PTree = PrincipalTree.PTree;
    edges = PrincipalTree.edges;
else
    projection = PrincipalTree.projection;
    PTree = PrincipalTree.PTree_post;
    edges = PrincipalTree.edges_post;
end

if ~isfield(params, 'FaceColor')
    Colors = [
        0.6350, 0.0780, 0.1840
        0.4660, 0.6740, 0.1880
        0, 0.4470, 0.7410
        0.8500, 0.3250, 0.0980
        0.9290, 0.6940, 0.1250
        0.3010, 0.7450, 0.9330
        0.4940, 0.1840, 0.5560];
    FaceColor = Colors(mod(1:length(levels),7),:);
else
    FaceColor = params.FaceColor;
end

if ~isfield(params, 'order')
    order = 1:length(levels);
else
    order = params.order;
end

annotations = order(annotations);
levels = levels(order);

%% plot the principal tree with the annotations
plotDDRtree(projection, PTree, edges, annotations, levels, FaceColor);

if params.prog_flag~=0
    ProgPath = PrincipalTree.ProgPath;
    [m, n] = size(ProgPath);
    ProgLevel = {'Others','Extracted path'};
    FaceColor = [0.4660, 0.6740, 0.1880
        0.8500, 0.3250, 0.0980];
    for i=1:m
        for j=1:n
            ProgL = ones(size(annotations));
            ProgL(ProgPath{i,j}) = 2;
            plotDDRtree(projection, PTree, edges, ProgL, ProgLevel, FaceColor);
        end
    end
end
end