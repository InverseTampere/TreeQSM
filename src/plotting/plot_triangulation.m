function plot_triangulation(QSM,fig,nf,AllTree)

% Plots the triangulation model of the stem's bottom part and the cylinder
% model (rest of the stem or the rest of the tree). The optional inputs
% "fig", "nf", "All" are the figure number, number of facets for the 
% cylinders, and if All = 1, then all the tree is plotted. 

n = nargin;
if n < 4
    AllTree = 0;
    if n < 3
        nf = 20;
        if n == 1
            fig = 1;
        end
    end
end

Vert = double(QSM.triangulation.vert);
Facets = double(QSM.triangulation.facet);
CylInd = QSM.triangulation.cylind;
fvd = QSM.triangulation.fvd;
Bran = QSM.cylinder.branch;
nc = size(Bran,1);
ind = (1:1:nc)';
C = ind(Bran == 1);
figure(fig)
patch('Vertices',Vert,'Faces',Facets,'FaceVertexCData',fvd,'FaceColor','flat')
hold on
if AllTree
    Ind = (CylInd:1:nc)';
else
    Ind = (CylInd:1:C(end))';
end
plot_cylinder_model(QSM.cylinder,fig,nf,1,'branch',Ind)
axis equal
hold off
alpha(1)
