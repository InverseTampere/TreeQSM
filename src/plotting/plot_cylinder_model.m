function plot_cylinder_model(cylinder,Color,fig,nf,alp,Ind)

% ---------------------------------------------------------------------
% PLOT_CYLINDER_MODEL.M       Plots the given cylinder model
%
% Version 1.2.0
% Latest update     3 Aug 2021
%
% Copyright (C) 2013-2021 Pasi Raumonen
% ---------------------------------------------------------------------

% Plots the cylinder model.
% cylinder  Structure array containin the cylinder info
%               (radius, length, start, axis, BranchOrder)
% fig       Figure number
% nf        Number of facets in the cyliders (in the thickest cylinder,
%               scales down with radius to 4 which is the minimum)
% alp       Alpha value (1 = no trasparency, 0 = complete transparency)
% Color     If equals to "order", colors the cylinders based on branching
%               order, otherwise colors each branch with unique color
% Ind       Indexes of cylinders to be plotted from a subset of cylinders
%               (Optional, if not given then all cylinders are plotted)

% Changes from version 1.1.0 to 1.2.0, 3 Aug 2021:
% 1) Changed the surface plot ("patch") so that the edges are not plotted
%    with separate color, so the surface looks more smooth. Similarly added
%    shading. (These are added at the end of the file)
% 2) Added cylinder branch "Bran" and branch order "BOrd" vectors where the
%    coloring options are defined to prevent some errors

% Changes from version 1.0.0 to 1.1.0, 13 July 2020:
% 1) Added option for choosing the coloring based either on branch order or
%    unique color for each branch
% 2) Removed the possibility of the input "cylinder" being a matrix
% 3) Added default values for inputs

n = nargin;
if n < 5
  alp = 1;
  if n < 4
    nf = 20;
    if n < 3
      fig = 1;
      if n == 1
        Color = 'order';
      end
    end
  end
end

Rad = cylinder.radius;
Len = cylinder.length;
Sta = cylinder.start;
%Sta = Sta-Sta(1,:);
Axe = cylinder.axis;
if strcmp(Color,'order')
  BOrd = cylinder.BranchOrder;
  Bran = cylinder.branch;
end
if strcmp(Color,'branch')
  Bran = cylinder.branch;
  BOrd = cylinder.BranchOrder;
end

if nargin == 6
  Rad = Rad(Ind);
  Len = Len(Ind);
  Sta = Sta(Ind,:);
  Axe = Axe(Ind,:);
  BOrd = BOrd(Ind);
  if strcmp(Color,'branch')
    Bran = Bran(Ind);
  end
end

nc = size(Rad,1); % Number of cylinder

if strcmp(Color,'order')
  Color = 1;
  % Color the cylinders in branches based on the branch order
  col = [
    0.00  0.00  1.00
    0.00  0.50  0.00
    1.00  0.00  0.00
    0.00  0.75  0.75
    0.75  0.00  0.75
    0.75  0.75  0.00
    0.25  0.25  0.25
    0.75  0.25  0.25
    0.95  0.95  0.00
    0.25  0.25  0.75
    0.75  0.75  0.75
    0.00  1.00  0.00
    0.76  0.57  0.17
    0.54  0.63  0.22
    0.34  0.57  0.92
    1.00  0.10  0.60
    0.88  0.75  0.73
    0.10  0.49  0.47
    0.66  0.34  0.65
    0.99  0.41  0.23];
  col = repmat(col,[10,1]);
elseif strcmp(Color,'branch')
  Color = 0;
  % Color the cylinders in branches with an unique color of each branch
  N = double(max(Bran));
  col = rand(N,3);
  Par = cylinder.parent;
  for i = 2:nc
    if Par(i) > 0 && Bran(Par(i)) ~= Bran(i)
      C = col(Bran(Par(i)),:);
      c = col(Bran(i),:);
      while sum(abs(C-c)) < 0.2
        c = rand(1,3);
      end
      col(Bran(i),:) = c;
    end
  end
end

Cir = cell(nf,2);
for i = 4:nf
  B = [cos((1/i:1/i:1)*2*pi)' sin((1/i:1/i:1)*2*pi)' zeros(i,1)];
  T = [cos((1/i:1/i:1)*2*pi)' sin((1/i:1/i:1)*2*pi)' ones(i,1)];
  Cir{i,1} = [B; T];
  Cir{i,2} = [(1:1:i)' (i+1:1:2*i)' [(i+2:1:2*i)'; i+1] [(2:1:i)'; 1]];
end

Vert = zeros(2*nc*(nf+1),3);
Facets = zeros(nc*(nf+1),4);
fvd = zeros(nc*(nf+1),3);
t = 1;
f = 1;

% Scale, rotate and translate the standard cylinders
for i = 1:nc
  n = ceil(sqrt(Rad(i)/Rad(1))*nf);
  n = min(n,nf);
  n = max(n,4);
  C = Cir{n,1};
  % Scale
  C(:,1:2) = Rad(i)*C(:,1:2);
  C(n+1:end,3) = Len(i)*C(n+1:end,3);
  % Rotate
  ang = real(acos(Axe(i,3)));
  Axis = cross([0 0 1]',Axe(i,:)');
  Rot = rotation_matrix(Axis,ang);
  C = (Rot*C')';
  % Translate
  C = mat_vec_subtraction(C,-Sta(i,:));
  Vert(t:t+2*n-1,:) = C;
  Facets(f:f+n-1,:) = Cir{n,2}+t-1;
  if Color == 1
    fvd(f:f+n-1,:) = repmat(col(BOrd(i)+1,:),[n 1]);
  else
    fvd(f:f+n-1,:) = repmat(col(Bran(i),:),[n 1]);
  end
  t = t+2*n;
  f = f+n;
end
t = t-1;
f = f-1;
Vert = Vert(1:t,:);
Facets = Facets(1:f,:);
fvd = fvd(1:f,:);

figure(fig)
plot3(Vert(1,1),Vert(1,2),Vert(1,3))
patch('Vertices',Vert,'Faces',Facets,'FaceVertexCData',fvd,'FaceColor','flat')
alpha(alp)
axis equal
grid on
view(-37.5,30)

shading flat
lightangle(gca,-45,30)
lighting gouraud