function [beta,sigma] = RegMixEM(x,y,k)
%REGMIXEM 2-D mixture of linear regression EM algorithm
%Reference: Benaglia T, Chauveau D, Hunter D, et al. mixtools: An r package for analyzing finite mixture models[J]. Journal of Statistical Software, 2009, 32(6): 1-29.
[nDim,nSamples] = size(x);
X = [x;ones(1,nSamples)];
beta = zeros(nDim+1,k);
sigma = ones(1,k);
mixCoefficients = ones(1,k) / k;
P = zeros(nSamples,k);
y = y';

%% Initialization
diff = 1;
t = 0;
maxiter = 40;
if nDim == 1 && k == 2;
    % initialize the beta with two simple lines
    [xmin,indmin] = min(x);
    ymin = y(indmin);
    [xmax,indmax] = max(x);
    ymax = y(indmax);
    xmedian = x(round(nSamples/2));
    ymedian = y(round(nSamples/2));
    beta(1,1) = (ymin - ymedian) / (xmin - xmedian);
    beta(2,1) = ymin - beta(1,1) * xmin;
    beta(1,2) = (ymax - ymedian) / (xmax - xmedian);
    beta(2,2) = ymax - beta(1,2) * xmax;
    % initialize the sigma by 2-order curve fitting
    [~,S] = polyfit(x,y',2);
    sigma = sigma * S.normr^2;
end
%% EM step
while(diff>1e-3)
    t = t + 1;
    %% Plot for 2D
    if nDim == 1 && t == (maxiter+1) 
        figure(t)
        plot(x,y','b*')
        hold on
        for j = 1:k
            xplot = min(x):0.01:max(x);
            plot(xplot,xplot*beta(1,j)+beta(2,j),'r-');
        end
        hold off
        xlim([min(x),max(x)])
        ylim([min(y),max(y)])
    end
    for i = 1:nSamples
        for j = 1:k
           P(i,j) = -(sum((y(i) - X(:,i)'*beta(:,j)).^2)) / sigma(j); %TO DO
        end
        P(i,:) = P(i,:) - max(P(i,:));
    end
    for i = 1:nSamples
       P(i,:) = P(i,:) + log(mixCoefficients); 
    end
    P = exp(P);
    for i = 1:nSamples
       P(i,:) = P(i,:) / sum(P(i,:)); 
    end
    for j = 1:k
        Wj = diag(P(:,j));
        beta(:,j) = inv(X*Wj*X')*X*Wj*y;
        sigma(j) = sum((sqrt(Wj)*(y - X'*beta(:,j))).^2) / trace(Wj);
    end
    for j = 1:k
        mixCoefficients(j) = sum(P(:,j));
    end
    mixCoefficients = mixCoefficients / sum(mixCoefficients);
    if t >= maxiter 
       break
    end
end


end

