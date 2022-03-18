function plotDDRtree(projection, PTree, edges, label, label_level, FaceColor)
%% ========================================================================
% Plot a principal tree
%
%--------------------------------------------------------------------------
% Input
%   projection  : projection of data points in the sapce of reduced dimension
%   PTree       : data points of the principal tree
%   edges       : edges of the principal tree
%   label       : annotations of data points
%   label_level : annotation levels
%   FaceColor   : facecolor of data points
%--------------------------------------------------------------------------
% Author: Lu Li
% update history: 08/10/2020
%% ========================================================================

%% initializations
Y = label;
U = unique(Y);
n_class = length(U);

%plot smaples
figure,
hold on

for i=1:n_class
    plot3(projection(1,Y==U(i)), projection(2,Y==U(i)), projection(3,Y==U(i)),...
        'o','MarkerFaceColor',FaceColor(i,:),'MarkerSize',10,'MarkerEdgeColor','k');
end
TreePoints = PTree(1:3,:);

% plot principal tree
[m,n] = size(edges);
for i=1:m
    for j=1:n
        if edges(i,j)~=0
            plot3([TreePoints(1,i), TreePoints(1,j)], [TreePoints(2,i), TreePoints(2,j)],...
                [TreePoints(3,i), TreePoints(3,j)],'-k','linewidth',4);
        end
    end
end
xlabel('DDR1')
ylabel('DDR2')
zlabel('DDR3')
legend(label_level);
grid;
set(gca,'FontSize',14);
view(-122,14)
end