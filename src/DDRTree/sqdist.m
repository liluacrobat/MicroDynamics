function d = sqdist(a,b)
% a, b : [D, N] = size(a)
% calculate the square distance between a, b

aa = sum(a.^2,1); 
bb = sum(b.^2,1); 
ab = a'*b;
d = abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
end
