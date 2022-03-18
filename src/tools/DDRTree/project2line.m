function [projection, theta] = project2line(x1, x2, x)
% Project data in x onto line defined by x1-->x2
% x1 : Feature-by-1 direction vector
% x2 : Feature-by-1 direction vector
% x : Feature-by-Sample data matrix


v = x2 - x1;
projection = zeros(size(x));
theta = zeros(1, size(x,2));
for n = 1:size(x,2)
    theta(n) = v'*(x(:,n)-x1)/(v'*v);
    projection(:,n) = v*theta(n) + x1;
end


