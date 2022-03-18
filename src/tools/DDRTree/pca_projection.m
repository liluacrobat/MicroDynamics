function [W,vs]= pca_projection(C, L)
% solve the problem size(C) = NxN, size(W) = NxL
% max_W trace( W' C W ) : W' W = I

[U,V] = eig(C);
[vs, eig_idx] = sort(diag(V),'descend');
W = U(:,eig_idx(1:L));