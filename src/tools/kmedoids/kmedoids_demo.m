close all; clear;
d = 2;
k = 3;
n = 500;
[X,label] = kmeansRnd(d,k,n);
y = kmedoids(X,k);
plotClass(X,label);
figure;
plotClass(X,y);