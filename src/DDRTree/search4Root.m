function root_index = search4Root(Pcurve, Data, subtype, normal_label)
%% ====================================================
% Search for root node on Pcurve
% =====================================================

MST = buildMST(Pcurve);
degree = sum(full(MST)>0, 1);
singleton = find(degree==1);
% find the 10 nearest samples of the leaf points
IDX = knnsearch(Data, Pcurve(singleton,:), 'K', 10);
% determine leaf points within HC sample
vote=zeros(size(singleton));
for i=1:length(vote)
    id = IDX(i,:);
    vote(i)=sum(strcmp(subtype(id), normal_label))/length(id);
end
root_index=singleton(vote>=0.7);
% determine the root point (winthin HC and have the largest distance to other leaf ponts)
Not_HC = setdiff(singleton,root_index);
root_dis=zeros(length(root_index),length(Not_HC));
for i=1:length(root_index)
    for j=1:length(Not_HC)
        [root_dis(i,j),  ~, ~] = graphshortestpath(MST, root_index(i), Not_HC(j));
    end
end
[~,idx]=sort(sum(root_dis,2),'descend');
root_index=root_index(idx);
end