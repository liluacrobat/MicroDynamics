function MST = buildMST(Pcurve)
% Build a MST from Pcurve
% Input:
%       Pcurve
% Output:
%       MST

G = sparse(squareform(pdist(Pcurve)));
[Tree, ~] = graphminspantree(G);
MST = Tree+Tree';
end
