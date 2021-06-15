function [CurveLength,Error_mse,rang,params_opt] = ElbowDDR(X,Para4ddr,f_sig,rang,num)
rang = sort(rang);
n = length(rang);
CurveLength = zeros(1, n);
Error_mse = zeros(1, n);
for i=1:length(rang)
    if f_sig==1
        Para4ddr.sigma=rang(i);
    else
        Para4ddr.lambda=rang(i);
    end
    [~, ~,~, ~, history] = DDRTree(X, Para4ddr);

    Error_mse(i)= history.mse(end);
    CurveLength(i) = history.length(end);
end
params_opt = ElbowPosition(Error_mse, CurveLength,rang,1, f_sig,10);
% save(savename,'Error_sorted_mse', 'CurveLength_sorted', 'rang_sorted');
end