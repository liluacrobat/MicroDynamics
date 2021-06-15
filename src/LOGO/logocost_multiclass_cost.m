function cost = logocost_multiclass_cost(alpha, P, Z)

%Function ComputeAlpha: compute alpha using linear search
Weight = P.Weight;
descent = P.descent;
lambda = P.lambda;

W = Weight-alpha*descent;
V = (W.^2)';
a = V*Z; %margin

index = isinf(exp(-a));
Result = log(1+exp(-a));
index = find(index==1);
Result(index) = -a(index);
cost = sum(Result) + lambda*sum(V);


%---------------------
return 

