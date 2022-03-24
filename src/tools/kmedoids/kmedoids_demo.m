close all; clear;
rng default
d = 2;
k = 3;
n = 500;

[X,label] = kmeansRnd(d,k,n);

[mapped_data,~,power]=compute_mapping(X','PCA',2);

y = kmedoids(X',k);
% plotClass(X,label);
figure;
plotClass(mapped_data',label);
figure;
plotClass(mapped_data',y);
1
