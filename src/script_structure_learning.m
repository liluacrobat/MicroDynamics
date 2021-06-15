function PrincipalTree = script_structure_learning(Feature_Table, Table_clinic, params)
%% ========================================================================
% Embedded structure learning
%
%--------------------------------------------------------------------------
% Input
%   Feature_Table : Table of the selected OTUs
%       -- tax  
%            Taxonomy of the selected OTUs
%       -- logRel 
%            Relative abundance of the selected OTUs after 10-base log 
%            Transformation
%       -- weight
%            Feature weight of the selected OTUs learned from LOGO
%   Table_clinic : Meta table of clinic information
% 
%   params       : Parameters
%       -- sigma
%            Bandwidth parameter
%       -- lambda
%            Regularization parameter for inverse graph embedding
%       -- col_label
%            The Column of the clinical information used for feature
%            selection
%--------------------------------------------------------------------------
% Output
%   PrincipalTree : Principal tree learned from data
%       -- W
%            Projection matrix
%       -- projection 
%            Projection of data points in the sapce of reduced dimension
%       -- edges
%            Edges of the principal tree
%       -- edges_post
%            Edges of the principal tree after post processing
%       -- PTree
%            Data points of the principal tree
%       -- PTree_post
%            Data points of the principal tree after post processing
%       -- ProgPath
%            Extracted progression path including the ordered index of samples
%       -- ProgPath_tree
%            Extracted progression path including the ordered index of tree
%            points
%       -- ProgDis
%            Progression distance of samples along the extracted path
%       -- DataProjection
%            Projection of data on the principalk tree
%--------------------------------------------------------------------------
% Author: Lu Li
% update history: 08/10/2020
%% ========================================================================

%% Initializations
para4ddr = params;
para4ddr.maxIter = 100; % maximum iterations
para4ddr.eps = 1e-9;    % relative objective difference
para4ddr.dim = 3;       % reduced dimension
para4ddr.gamma = 2;     % regularization parameter for k-means

%% Learn the principal tree
[W, projection, edges, PTree, ~] = DDRTree(Feature_Table.logRel, para4ddr);

PrincipalTree.W = W;
PrincipalTree.projection  = projection;
PrincipalTree.edges  = edges;
PrincipalTree.PTree  = PTree;

subtype = table2cell(Table_clinic(:,params.col_label));
[extracted_path, extracted_pathDist, PTree_etd, extracted_curve,...
    edges_post, DataProjection] = StruPostProcess(PTree, projection, subtype);

PrincipalTree.PTree_post = PTree_etd;
PrincipalTree.ProgPath = extracted_path;
PrincipalTree.ProgPath_tree = extracted_curve;
PrincipalTree.ProgDis = extracted_pathDist;
PrincipalTree.edges_post = edges_post;
PrincipalTree.DataProjection = DataProjection;
end