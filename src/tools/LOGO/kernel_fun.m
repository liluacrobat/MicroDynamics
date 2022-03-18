function prob = kernel_fun(dist,Para)
[~,I] = sort(dist);
prob = zeros(size(I));
kernel = Para.kernel;
dis_neighbor = dist(I(1:Para.sigma));
dis_neighbor = dis_neighbor/max(dis_neighbor);
switch kernel
    case 'triangular'
        prob(I(1:Para.sigma)) = 1-abs(dis_neighbor);
    case 'parabolic'
        prob(I(1:Para.sigma)) = 1-(dis_neighbor).^2;
    case 'quartic'
        prob(I(1:Para.sigma)) = (1-(dis_neighbor).^2).^2;
    case 'triweight'
        prob(I(1:Para.sigma)) = (1-(dis_neighbor).^2).^3;
    case 'tricube'
        prob(I(1:Para.sigma)) = (1-(abs(dis_neighbor)).^3).^3;
    case 'cosine'
        prob(I(1:Para.sigma)) = cos(pi/2*dis_neighbor);
    otherwise
        prob(I(1:Para.sigma)) = 1;
end
prob = prob/sum(prob);
%% ==================End of the code===================================
end
