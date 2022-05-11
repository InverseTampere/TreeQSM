% This file is part of TREEQSM.
%
% TREEQSM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% TREEQSM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with TREEQSM.  If not, see <http://www.gnu.org/licenses/>.

function triangulation = curve_based_triangulation(P,TriaHeight,TriaWidth)

% ---------------------------------------------------------------------
% CURVE_BASED_TRIANGULATION.M   Reconstructs a triangulation for the
%                               stem-buttress surface based on boundary curves
%
% Version 1.1.0
% Latest update     3 May 2022
%
% Copyright (C) 2015-2022 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Inputs:
% P             Point cloud of the stem to be triangulated
% TriaHeight    Height of the triangles
% TriaWidth     Width of the triangles
%
% Output:
% triangulation  Structure field defining the triangulation. Contains
%                   the following main fields:
%   vert            Vertices of the triangulation model (nv x 3)-matrix
%   facet           Facets (triangles) of the triangulation 
%                     (the vertices forming the facets)
%   fvd             Color information of the facets for plotting with "patch"
%   volume          Volume enclosed by the facets in liters
%   bottom          The z-coordinate of the bottom of the model
%   top             The z-coordinate of the top of the model
%   triah           TriaHeight
%   triaw           TriaWidth
% ---------------------------------------------------------------------

% Changes from version 1.0.2 to 1.1.0, 3 May 2022:
% 1) Increased the radius of the balls at seed points from TriaWidth to 
%    2*TriaWidth in the input of "boundary_curve"
% 2) Added triangle orientation check after the side is covered with
%    triangles so that the surface normals are pointing outward 
% 3) Modified the check if the new boundary curve changes only a little and 
%    then stop reconstruction  
% 4) Added halving the triangle height if the boundary curve length has
%    increased three times.
% 5) Changed the bottom level from the smallest z-coordinate to the  
%    average of the lowest 100 z-coordinates. 
% 6) Minor streamlining the code and added more comments

% Changes from version 1.0.2 to 1.0.3, 11 Aug 2020:
% 1) Small changes in the code when computing the delaunay triangulation
%    of the top layer

% Changes from version 1.0.1 to 1.0.2, 15 Jan 2020:
% 1) Added side surface areas (side, top, bottom) to output as fields

% Changes from version 1.0.0 to 1.0.1, 26 Nov 2019:
% 1) Removed the plotting of the triangulation model at the end of the code

%% Determine the first boundary curve
np = size(P,1);
[~,I] = sort(P(:,3),'descend');
P = P(I,:);
Hbot = mean(P(end-100:end,3));
Htop = P(1,3);
N = ceil((Htop-Hbot)/TriaHeight);
Vert = zeros(1e5,3);
Tria = zeros(1e5,3);
TriaLay = zeros(1e5,1);
VertLay = zeros(1e5,1,'uint16');
Curve = zeros(0,3);
i = 0; % the layer whose cross section is under reconstruction
ps = 1;
while P(ps,3) > Htop-i*TriaHeight
  ps = ps+1;
end
pe = ps;
while i < N/4 && isempty(Curve)
  % Define thin horizontal cross section of the stem
  i = i+1;
  ps = pe+1;
  k = 1;
  while P(ps+k,3) > Htop-i*TriaHeight
    k = k+1;
  end
  pe = ps+k-1;
  PSection = P(ps:pe,:);

  % Create initial boundary curve:
  iter = 0;
  while iter <= 15 && isempty(Curve)
    iter = iter+1;
    Curve = initial_boundary_curve(PSection,TriaWidth);
  end
end

if isempty(Curve)
  triangulation = zeros(0,1);
  disp('  No triangulation: Problem with the first curve')
  return
end

% make the height of the curve even:
Curve(:,3) = max(Curve(:,3));
% Save vertices:
nv = size(Curve,1); % number of vertices in the curve
Vert(1:nv,:) = Curve;
VertLay(1:nv) = i;
t = 0;
m00 = size(Curve,1);

%% Determine the other boundary curves and the triangulation downwards
i0 = i;
i = i0+1;
nv0 = 0;
LayerBottom = Htop-i*TriaHeight;
while i <= N && pe < np
  %% Define thin horizontal cross section of the stem
  ps = pe+1;
  k = 1;
  while ps+k <= np && P(ps+k,3) > LayerBottom
    k = k+1;
  end
  pe = ps+k-1;
  PSection = P(ps:pe,:);

  %% Create boundary curves using the previous curves as seeds
  if i > i0+1
    nv0 = nv1;
  end
  % Define seed points:
  Curve(:,3) = Curve(:,3)-TriaHeight;
  Curve0 = Curve;

  % Create new boundary curve
  [Curve,Ind] = boundary_curve(PSection,Curve,2*TriaWidth,1.5*TriaWidth);

  if isempty(Curve)
    disp('  No triangulation: Empty curve')
    triangulation = zeros(0,1);
    return
  end
  Curve(:,3) = max(Curve(:,3));

  %% Check if the curve intersects itself
  [Intersect,IntersectLines] = check_self_intersection(Curve(:,1:2));

  %% If self-intersection, try to modify the curve
  j = 1;
  while Intersect && j <= 10
    n = size(Curve,1);
    CrossLines = (1:1:n)';
    NumberOfIntersections = cellfun('length',IntersectLines(:,1));
    I = NumberOfIntersections > 0;
    CrossLines = CrossLines(I);
    CrossLen = vertcat(IntersectLines{I,2});
    if length(CrossLen) == length(CrossLines)
      LineEle = Curve([2:end 1],:)-Curve(1:end,:);
      d = sqrt(sum(LineEle.*LineEle,2));
      m = length(CrossLines);
      for k = 1:2:m
        if CrossLines(k) ~= n
          Curve(CrossLines(k)+1,:) = Curve(CrossLines(k),:)+...
            0.9*CrossLen(k)/d(CrossLines(k))*LineEle(CrossLines(k),:);
        else
          Curve(1,:) = Curve(CrossLines(k),:)+...
            0.9*CrossLen(k)/d(CrossLines(k))*LineEle(CrossLines(k),:);
        end
      end
      [Intersect,IntersectLines] = check_self_intersection(Curve(:,1:2));
      j = j+1;
    else
      j = 11;
    end
  end

  m = size(Curve,1);
  if Intersect
    %% Curve self-intersects, use previous curve to extrapolate to the bottom
    H = Curve0(1,3)-Hbot;
    if H > 0.75 && Intersect
      triangulation = zeros(0,1);
      disp(['  No triangulation: Self-intersection at ',...
        num2str(H),' m from the bottom'])
      return
    end
    Curve = Curve0;
    Curve(:,3) = Curve(:,3)-TriaHeight;
    Nadd = floor(H/TriaHeight)+1;
    m = size(Curve,1);
    Ind = [(1:1:m)' [(2:1:m)'; 1]];
    T = H/Nadd;
    for k = 1:Nadd
      if k > 1
        Curve(:,3) = Curve(:,3)-T;
      end
      Vert(nv+1:nv+m,:) = Curve;
      VertLay(nv+1:nv+m) = i;
      %% Define the triangulation between two boundary curves
      nv1 = nv;
      nv = nv+m;
      t0 = t+1;
      pass = false;
      for j = 1:m
        if Ind(j,2) > 0 && j < m
          t = t+1;
          Tria(t,:) = [nv1+j nv0+Ind(j,:)];
          t = t+1;
          Tria(t,:) = [nv1+j nv0+Ind(j,2) nv1+j+1];
        elseif Ind(j,2) > 0 && ~pass
          t = t+1;
          Tria(t,:) = [nv1+j nv0+Ind(j,:)];
          t = t+1;
          Tria(t,:) = [nv1+j nv0+Ind(j,2) nv1+1];
        elseif Ind(j,2) == 0 && j < m
          t = t+1;
          Tria(t,:) = [nv1+j nv0+Ind(j,1) nv1+j+1];
        elseif Ind(j,2) == 0 && ~pass
          t = t+1;
          Tria(t,:) = [nv1+j nv0+Ind(j,1) nv1+1];
        elseif j == 1 && Ind(j,2) == -1
          t = t+1;
          Tria(t,:) = [nv nv1 nv0+1];
          t = t+1;
          Tria(t,:) = [nv nv0+1 nv1+1];
          t = t+1;
          Tria(t,:) = [nv0+1 nv0+2 nv1+1];
          t = t+1;
          Tria(t,:) = [nv1+1 nv0+2 nv0+3];
          t = t+1;
          Tria(t,:) = [nv1+1 nv0+3 nv1+2];
          pass = true;
        elseif Ind(j,2) == -1 && j < m
          t = t+1;
          Tria(t,:) = [nv1+j nv0+Ind(j,1) nv0+Ind(j,1)+1];
          t = t+1;
          Tria(t,:) = [nv1+j nv0+Ind(j,1)+1 nv1+j+1];
          t = t+1;
          Tria(t,:) = [nv0+Ind(j,1)+1 nv0+Ind(j,1)+2 nv1+j+1];
        elseif Ind(j,2) == -1 && ~pass
          t = t+1;
          Tria(t,:) = [nv1+j nv0+Ind(j,1) nv0+Ind(j,1)+1];
          t = t+1;
          Tria(t,:) = [nv1+j nv0+Ind(j,1)+1 nv1+1];
          t = t+1;
          Tria(t,:) = [nv0+Ind(j,1)+1 nv0+1 nv1+1];
        end
      end

      TriaLay(t0:t) = i;
      i = i+1;
      nv0 = nv1;
    end
    i = N+1;

  else
    %% No self-intersection, proceed with triangulation and new curves
    Vert(nv+1:nv+m,:) = Curve;
    VertLay(nv+1:nv+m) = i;

    %% If little change between Curve and Curve0, stop the reconstruction
    C = intersect(Curve0,Curve,"rows");
    if size(C,1) > 0.7*size(Curve,1)
      N = i;
    end

    %% If the boundary curve has grown much longer than originally, then
    % decrease the triangle height
    if m > 3*m00
      TriaHeight = TriaHeight/2; % use half the height
      N = N+ceil((N-i)/2); % update the number of layers 
      m00 = m;
    end

    %% Define the triangulation between two boundary curves
    nv1 = nv;
    nv = nv+m;
    t0 = t+1;
    pass = false;
    for j = 1:m
      if Ind(j,2) > 0 && j < m
        t = t+1;
        Tria(t,:) = [nv1+j nv0+Ind(j,:)];
        t = t+1;
        Tria(t,:) = [nv1+j nv0+Ind(j,2) nv1+j+1];
      elseif Ind(j,2) > 0 && ~pass
        t = t+1;
        Tria(t,:) = [nv1+j nv0+Ind(j,:)];
        t = t+1;
        Tria(t,:) = [nv1+j nv0+Ind(j,2) nv1+1];
      elseif Ind(j,2) == 0 && j < m
        t = t+1;
        Tria(t,:) = [nv1+j nv0+Ind(j,1) nv1+j+1];
      elseif Ind(j,2) == 0 && ~pass
        t = t+1;
        Tria(t,:) = [nv1+j nv0+Ind(j,1) nv1+1];
      elseif j == 1 && Ind(j,2) == -1
        t = t+1;
        Tria(t,:) = [nv nv1 nv0+1];
        t = t+1;
        Tria(t,:) = [nv nv0+1 nv1+1];
        t = t+1;
        Tria(t,:) = [nv0+1 nv0+2 nv1+1];
        t = t+1;
        Tria(t,:) = [nv1+1 nv0+2 nv0+3];
        t = t+1;
        Tria(t,:) = [nv1+1 nv0+3 nv1+2];
        pass = true;
      elseif Ind(j,2) == -1 && j < m
        t = t+1;
        Tria(t,:) = [nv1+j nv0+Ind(j,1) nv0+Ind(j,1)+1];
        t = t+1;
        Tria(t,:) = [nv1+j nv0+Ind(j,1)+1 nv1+j+1];
        t = t+1;
        Tria(t,:) = [nv0+Ind(j,1)+1 nv0+Ind(j,1)+2 nv1+j+1];
      elseif Ind(j,2) == -1 && ~pass
        t = t+1;
        Tria(t,:) = [nv1+j nv0+Ind(j,1) nv0+Ind(j,1)+1];
        t = t+1;
        Tria(t,:) = [nv1+j nv0+Ind(j,1)+1 nv1+1];
        t = t+1;
        Tria(t,:) = [nv0+Ind(j,1)+1 nv0+1 nv1+1];
      end
    end

    TriaLay(t0:t) = i;
    i = i+1;
    LayerBottom = LayerBottom-TriaHeight;
  end

end
Vert = Vert(1:nv,:);
VertLay = VertLay(1:nv);
Tria = Tria(1:t,:);
TriaLay = TriaLay(1:t);

%% Check the orientation of the triangles 
% so that surface normals are outward pointing
a = round(t/10); % select the top triangles
U = Vert(Tria(1:a,2),:)-Vert(Tria(1:a,1),:);
V = Vert(Tria(1:a,3),:)-Vert(Tria(1:a,1),:);
Center = mean(Vert(1:nv-1,:)); % the center of the stem
C = Vert(Tria(1:a,1),:)+0.25*V+0.25*U;
W = C(:,1:2)-Center(1:2); % vectors from the triagles to the stem's center
Normals = cross(U,V);
if nnz(sum(Normals(:,1:2).*W,2) < 0) > 0.5*length(C)
  Tria(1:t,1:2) = [Tria(1:t,2) Tria(1:t,1)];
end

% U = Vert(Tria(1:t,2),:)-Vert(Tria(1:t,1),:);
% V = Vert(Tria(1:t,3),:)-Vert(Tria(1:t,1),:);
% Normals = cross(U,V);
% Normals = normalize(Normals);
% C = Vert(Tria(1:t,1),:)+0.25*V+0.25*U;
% fvd = ones(t,1);
% figure(5)
% point_cloud_plotting(P(1,:),5,6)
% patch('Vertices',Vert,'Faces',Tria,'FaceVertexCData',fvd,'FaceColor','flat')
% alpha(1)
% hold on
% arrow_plot(C,0.1*Normals,5)
% hold off
% axis equal
% pause


%% Remove possible double triangles
nt = size(Tria,1);
Keep = true(nt,1);
Scoord = Vert(Tria(:,1),:)+Vert(Tria(:,2),:)+Vert(Tria(:,3),:);
S = sum(Scoord,2);
[part,CC] = cubical_partition(Scoord,2*TriaWidth);
for j = 1:nt-1
  if Keep(j)
    points = part(CC(j,1)-1:CC(j,1)+1,CC(j,2)-1:CC(j,2)+1,CC(j,3)-1:CC(j,3)+1);
    points = vertcat(points{:});
    I = S(j) == S(points);
    J = points ~= j;
    I = I&J&Keep(points);
    if any(I)
      p = points(I);
      I = intersect(Tria(j,:),Tria(p,:));
      if length(I) == 3
        Keep(p) = false;
      end
    end
  end
end
Tria = Tria(Keep,:);
TriaLay = TriaLay(Keep);


%% Generate triangles for the horizontal layers and compute the volumes
% Triangles of the ground layer
% Select the boundary curve:
N = double(max(VertLay));
I = VertLay == N;
Vert(I,3) = Hbot;
ind = (1:1:nv)';
ind = ind(I);
Curve = Vert(I,:); % Boundary curve of the bottom
n = size(Curve,1);
if n < 10
  triangulation = zeros(0,1);
  disp('  No triangulation: Ground layer boundary curve too small')
  return
end

% Define Delaunay triangulation for the bottom
C = zeros(n,2);
C(:,1) = (1:1:n)';
C(1:n-1,2) = (2:1:n)';
C(n,2) = 1;
warning off
dt = delaunayTriangulation(Curve(:,1),Curve(:,2),C);
In = dt.isInterior();
GroundTria = dt(In,:);
Points = dt.Points;
warning on
if size(Points,1) > size(Curve,1)
  disp('  No triangulation: Problem with delaunay in the bottom layer')
  triangulation = zeros(0,1);
  return
end
GroundTria0 = GroundTria;
GroundTria(:,1) = ind(GroundTria(:,1));
GroundTria(:,2) = ind(GroundTria(:,2));
GroundTria(:,3) = ind(GroundTria(:,3));

% Compute the normals and areas
U = Curve(GroundTria0(:,2),:)-Curve(GroundTria0(:,1),:);
V = Curve(GroundTria0(:,3),:)-Curve(GroundTria0(:,1),:);
Cg = Curve(GroundTria0(:,1),:)+0.25*V+0.25*U;
Ng = cross(U,V);
I = Ng(:,3) > 0; % Check orientation
Ng(I,:) = -Ng(I,:);
Ag = 0.5*sqrt(sum(Ng.*Ng,2));
Ng = 0.5*[Ng(:,1)./Ag Ng(:,2)./Ag Ng(:,3)./Ag];

% Remove possible negative area triangles:
I = Ag > 0;   Ag = Ag(I);   Cg = Cg(I,:);   Ng = Ng(I,:);
GroundTria = GroundTria(I,:);

% Update the triangles:
Tria = [Tria; GroundTria];
TriaLay = [TriaLay; (N+1)*ones(size(GroundTria,1),1)];

if abs(sum(Ag)-polyarea(Curve(:,1),Curve(:,2))) > 0.001*sum(Ag)
  disp('  No triangulation: Problem with delaunay in the bottom layer')
  triangulation = zeros(0,1);
  return
end

% Triangles of the top layer
% Select the top curve:
N = double(min(VertLay));
I = VertLay == N;
ind = (1:1:nv)';
ind = ind(I);
Curve = Vert(I,:);
CenterTop = mean(Curve);
%  Delaunay triangulation of the top:
n = size(Curve,1);
C = zeros(n,2);
C(:,1) = (1:1:n)';
C(1:n-1,2) = (2:1:n)';
C(n,2) = 1;
warning off
dt = delaunayTriangulation(Curve(:,1),Curve(:,2),C);
Points = dt.Points;
warning on
if min(size(dt)) == 0 || size(Points,1) > size(Curve,1)
  disp('  No triangulation: Problem with delaunay in the top layer')
  triangulation = zeros(0,1);
  return
end
In = dt.isInterior();
TopTria = dt(In,:);
TopTria0 = TopTria;
TopTria(:,1) = ind(TopTria(:,1));
TopTria(:,2) = ind(TopTria(:,2));
TopTria(:,3) = ind(TopTria(:,3));

% Compute the normals and areas:
U = Curve(TopTria0(:,2),:)-Curve(TopTria0(:,1),:);
V = Curve(TopTria0(:,3),:)-Curve(TopTria0(:,1),:);
Ct = Curve(TopTria0(:,1),:)+0.25*V+0.25*U;
Nt = cross(U,V);
I = Nt(:,3) < 0;
Nt(I,:) = -Nt(I,:);
At = 0.5*sqrt(sum(Nt.*Nt,2));
Nt = 0.5*[Nt(:,1)./At Nt(:,2)./At Nt(:,3)./At];

% Remove possible negative area triangles:
I = At > 0;   At = At(I);   Ct = Ct(I,:);   Nt = Nt(I,:);
TopTria = TopTria(I,:);

% Update the triangles:
Tria = [Tria; TopTria];
TriaLay = [TriaLay; N*ones(size(TopTria,1),1)];

if abs(sum(At)-polyarea(Curve(:,1),Curve(:,2))) > 0.001*sum(At)
  disp('  No triangulation: Problem with delaunay in the top layer')
  triangulation = zeros(0,1);
  return
end

% Triangles of the side
B = TriaLay <= max(VertLay) & TriaLay > 1;
U = Vert(Tria(B,2),:)-Vert(Tria(B,1),:);
V = Vert(Tria(B,3),:)-Vert(Tria(B,1),:);
Cs = Vert(Tria(B,1),:)+0.25*V+0.25*U;
Ns = cross(U,V);
As = 0.5*sqrt(sum(Ns.*Ns,2));
Ns = 0.5*[Ns(:,1)./As Ns(:,2)./As Ns(:,3)./As];
I = As > 0;  Ns = Ns(I,:); As = As(I); Cs = Cs(I,:);

% Volumes in liters
VTotal = sum(At.*sum(Ct.*Nt,2))+sum(As.*sum(Cs.*Ns,2))+sum(Ag.*sum(Cg.*Ng,2));
VTotal = round(10000*VTotal/3)/10;

if VTotal < 0
  disp('  No triangulation: Problem with volume')
  triangulation = zeros(0,1);
  return
end

V = Vert(Tria(:,1),1:2)-CenterTop(1:2);
fvd = sqrt(sum(V.*V,2));
triangulation.vert = single(Vert);
triangulation.facet = uint16(Tria);
triangulation.fvd = single(fvd);
triangulation.volume = VTotal;
triangulation.SideArea = sum(As);
triangulation.BottomArea = sum(Ag);
triangulation.TopArea = sum(At);
triangulation.bottom = min(Vert(:,3));
triangulation.top = max(Vert(:,3));
triangulation.triah = TriaHeight;
triangulation.triaw = TriaWidth;

% figure(5)
% point_cloud_plotting(P,5,6)
% patch('Vertices',Vert,'Faces',Tria,'FaceVertexCData',fvd,'FaceColor','flat')
% % hold on
% % arrow_plot(Cs,0.2*Ns,5)
% % hold off
% % axis equal
% alpha(1)

