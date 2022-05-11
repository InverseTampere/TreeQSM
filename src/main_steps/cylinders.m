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

function cylinder = cylinders(P,cover,segment,inputs)

% ---------------------------------------------------------------------
% CYLINDERS.M       Fits cylinders to the branch-segments of the point cloud
%
% Version 3.0.0
% Latest update     1 Now 2018
%
% Copyright (C) 2013-2018 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Reconstructs the surface and volume of branches of input tree with
% cylinders. Subdivides each segment to smaller regions to which cylinders
% are fitted in least squares sense. Returns the cylinder information and
% in addition the child-relation of the cylinders plus the cylinders in
% each segment.
% ---------------------------------------------------------------------
% Inputs:
% P         Point cloud, matrix
% cover     Cover sets
% segment   Segments
% input     Input parameters of the reconstruction:
%   MinCylRad   Minimum cylinder radius, used in the taper corrections
%   ParentCor   Radius correction based on radius of the parent: radii in
%                 a branch are usually smaller than the radius of the parent
%                 cylinder in the parent branch
%   TaperCor    Parabola taper correction of radii inside branches.
%   GrowthVolCor  If 1, use growth volume correction
%   GrowthVolFac  Growth volume correction factor
%
% Outputs:
% cylinder  Structure array containing the following cylinder info:
%   radius        Radii of the cylinders, vector
%   length        Lengths of the cylinders, vector
%   axis          Axes of the cylinders, matrix
%   start         Starting points of the cylinders, matrix
%   parent        Parents of the cylinders, vector
%   extension     Extensions of the cylinders, vector
%   branch        Branch of the cylinder
%   BranchOrder   Branching order of the cylinder
%   PositionInBranch    Position of the cylinder inside the branch
%   mad           Mean absolute distances of points from the cylinder
%                           surface, vector
%   SurfCov       Surface coverage measure, vector
%   added         Added cylinders, logical vector
%   UnModRadius   Unmodified radii
% ---------------------------------------------------------------------

% Changes from version 3.0.0 to 3.1.0, 6 Oct 2021:
% 1) Added the growth volume correction option ("growth_volume_correction")
%    back, which was removed from the previous version by a mistake. The
%    "growth_volume_correction" function was also corrected.
% 2) Added the fields "branch", "BranchOrder", "PositionInBranch" to the
%    output structure "cylinder"
% 3) Removed the fields "CylsInSegment" and "ChildCyls" from the output
%    structure "cylinder"

% Changes from version 2.0.0 to 3.0.0, 13 Aug 2020:
% Many comprehensive and small changes:
% 1) "regions" and "cylinder_fitting" are combined into "cylinder_fitting"
%   and the process is more adaptive as it now fits at least 3 (up to 10)
%   cylinders of different lengths for each region.
% 2) "lcyl" and "FilRad" parameters are not used anymore
% 3) Surface coverage ("SurfCov") and mean absolute distance ("mad") are
%   added to the cylinder structure as fields.
% 4) Surface coverage filtering is used in the definition of the regions
%   and removing outliers
% 5) "adjustments" has many changes, particularly in the taper corrections
%   where the parabola-taper curve is fitted to all the data with surface
%   coverage as a weight. Adjustment of radii based on the parabola is
%   closer the parabola the smaller the surface coverage. For the stem the
%   taper correction is the same as for the branches. The minimum and
%   maximum radii corrections are also modified.
% 6) Syntax has changed, particularly for the "cyl"-structure

% Changes from version 2.1.0 to 2.1.1, 26 Nov 2019:
% 1) Increased the minimum number "n" of estimated cylinders for
%    initialization of vectors at the beginning of the code. This is done
%    to make sure that trees without branches will not cause errors.

% Changes from version 2.0.0 to 2.1.0, 3 Oct 2019:
% 1) Bug fix: UnmodRadius is now defined as it should, as the radius after
%    least squares fitting but without parent, taper or growth vol. corrections
% 2) Bug fix: Correction in "least_squares_cylinder.m", calculates the
%    starting point of the cylinder now correctly.
% 3) Bug fix: Correct errors related to combining data when a fitted
%    cylinder is replaced with two shorter ones, in "cylinder_fitting"
% 4) Removed some unnecessary command lines for computing radius estimates
%    in "regions"

%% Initialization of variables
Segs = segment.segments;
SPar = segment.ParentSegment;
SChi = segment.ChildSegment;
NumOfSeg = max(size(Segs));   % number of segments
n = max(2000,min(40*NumOfSeg,2e5));
c = 1;  % number of cylinders determined
CChi = cell(n,1); % Children of the cylinders
CiS = cell(NumOfSeg,1); % Cylinders in the segment
cylinder.radius = zeros(n,1,'single');
cylinder.length = zeros(n,1,'single');
cylinder.start = zeros(n,3,'single');
cylinder.axis = zeros(n,3,'single');
cylinder.parent = zeros(n,1,'uint32');
cylinder.extension = zeros(n,1,'uint32');
cylinder.added = false(n,1);
cylinder.UnmodRadius = zeros(n,1,'single');
cylinder.branch = zeros(n,1,'uint16');
cylinder.SurfCov = zeros(n,1,'single');
cylinder.mad = zeros(n,1,'single');

%% Determine suitable order of segments (from trunk to the "youngest" child)
bases = (1:1:NumOfSeg)';
bases = bases(SPar(:,1) == 0);
nb = length(bases);
SegmentIndex = zeros(NumOfSeg,1);
nc = 0;
for i = 1:nb
  nc = nc+1;
  SegmentIndex(nc) = bases(i);
  S = vertcat(SChi{bases(i)});
  while ~isempty(S)
    n = length(S);
    SegmentIndex(nc+1:nc+n) = S;
    nc = nc+n;
    S = vertcat(SChi{S});
  end
end

%% Fit cylinders individually for each segment
for k = 1:NumOfSeg
  si = SegmentIndex(k);
  if si > 0
    %% Some initialization about the segment
    Seg = Segs{si};      % the current segment under analysis
    nl = max(size(Seg));  % number of cover set layers in the segment
    [Sets,IndSets] = verticalcat(Seg); % the cover sets in the segment

    ns = length(Sets);   % number of cover sets in the current segment
    Points = vertcat(cover.ball{Sets}); % the points in the segments
    np = length(Points);         % number of points in the segment

    % Determine indexes of points for faster definition of regions
    BallSize = cellfun('length',cover.ball(Sets));
    IndPoints = ones(nl,2); % indexes for points in each layer of the segment
    for j = 1:nl
      IndPoints(j,2) = sum(BallSize(IndSets(j,1):IndSets(j,2)));
    end
    IndPoints(:,2) = cumsum(IndPoints(:,2));
    IndPoints(2:end,1) = IndPoints(2:end,1)+IndPoints(1:end-1,2);
    Base = Seg{1};          % the base of the segment
    nb = IndPoints(1,2); % number of points in the base

    % Reconstruct only large enough segments
    if nl > 1 && np > nb && ns > 2 && np > 20 && ~isempty(Base)

      %% Cylinder fitting
      [cyl,Reg] = cylinder_fitting(P,Points,IndPoints,nl,si);
      nc = numel(cyl.radius);

      %% Search possible parent cylinder
      if nc > 0 && si > 1
        [PC,cyl,added] = parent_cylinder(SPar,SChi,CiS,cylinder,cyl,si);
        nc = numel(cyl.radius);
      elseif si == 1
        PC = zeros(0,1);
        added = false;
      else
        added = false;
      end
      cyl.radius0 = cyl.radius;

      %% Modify cylinders
      if nc > 0
        % Define parent cylinder:
        parcyl.radius = cylinder.radius(PC);
        parcyl.length = cylinder.length(PC);
        parcyl.start = cylinder.start(PC,:);
        parcyl.axis = cylinder.axis(PC,:);
        % Modify the cylinders
        cyl = adjustments(cyl,parcyl,inputs,Reg);
      end

      %% Save the cylinders
      % if at least one acceptable cylinder, then save them
      Accept = nc > 0 & min(cyl.radius(1:nc)) > 0;
      if Accept
        % If the parent cylinder exists, set the parent-child relations
        if ~isempty(PC)
          cylinder.parent(c) = PC;
          if cylinder.extension(PC) == c
            I = cylinder.branch(PC);
            cylinder.branch(c:c+nc-1) = I;
            CiS{I} = [CiS{I}; linspace(c,c+nc-1,nc)'];
          else
            CChi{PC} = [CChi{PC}; c];
            cylinder.branch(c:c+nc-1) = si;
            CiS{si} = linspace(c,c+nc-1,nc)';
          end
        else
          cylinder.branch(c:c+nc-1) = si;
          CiS{si} = linspace(c,c+nc-1,nc)';
        end

        cylinder.radius(c:c+nc-1) = cyl.radius(1:nc);
        cylinder.length(c:c+nc-1) = cyl.length(1:nc);
        cylinder.axis(c:c+nc-1,:) = cyl.axis(1:nc,:);
        cylinder.start(c:c+nc-1,:) = cyl.start(1:nc,:);
        cylinder.parent(c+1:c+nc-1) = linspace(c,c+nc-2,nc-1);
        cylinder.extension(c:c+nc-2) = linspace(c+1,c+nc-1,nc-1);
        cylinder.UnmodRadius(c:c+nc-1) = cyl.radius0(1:nc);
        cylinder.SurfCov(c:c+nc-1) = cyl.SurfCov(1:nc);
        cylinder.mad(c:c+nc-1) = cyl.mad(1:nc);
        if added
          cylinder.added(c) = true;
          cylinder.added(c) = true;
        end
        c = c+nc; % number of cylinders so far (plus one)

      end
    end
  end
end
c = c-1; % number of cylinders


%% Define outputs
names = fieldnames(cylinder);
n = max(size(names));
for k = 1:n
  cylinder.(names{k}) = single(cylinder.(names{k})(1:c,:));
end
if c <= 2^16
  cylinder.parent = uint16(cylinder.parent);
  cylinder.extension = uint16(cylinder.extension);
end
nb = max(cylinder.branch);
if nb <= 2^8
  cylinder.branch = uint8(cylinder.branch);
elseif nb <= 2^16
  cylinder.branch = uint16(cylinder.branch);
end
cylinder.added = logical(cylinder.added);

% Define the branching order:
BOrd = zeros(c,1);
for i = 1:c
  if cylinder.parent(i) > 0
    p = cylinder.parent(i);
    if cylinder.extension(p) == i
      BOrd(i) = BOrd(p);
    else
      BOrd(i) = BOrd(p)+1;
    end
  end
end
cylinder.BranchOrder = uint8(BOrd);
% Define the cylinder position inside the branch
PiB = ones(c,1);
for i = 1:NumOfSeg
  C = CiS{i};
  if ~isempty(C)
    n = length(C);
    PiB(C) = (1:1:n)';
  end
end
if max(PiB) <= 2^8
  cylinder.PositionInBranch = uint8(PiB);
else
  cylinder.PositionInBranch = uint16(PiB);
end

% Growth volume correction
if inputs.GrowthVolCor && c > 0
  cylinder = growth_volume_correction(cylinder,inputs);
end

end % End of main function


function [cyl,Reg] = cylinder_fitting(P,Points,Ind,nl,si)

if nl > 6
  i0 = 1;     i = 4; % indexes of the first and last layers of the region
  t = 0;
  Reg = cell(nl,1);
  cyls = cell(11,1);
  regs = cell(11,1);
  data = zeros(11,4);
  while i0 < nl-2
    %% Fit at least three cylinders of different lengths
    bot = Points(Ind(i0,1):Ind(i0+1,2));
    Bot = average(P(bot,:)); % Bottom axis point of the region
    again = true;
    j = 0;
    while i+j <= nl && j <= 10 && (j <= 2 || again)
      %% Select points and estimate axis
      RegC = Points(Ind(i0,1):Ind(i+j,2)); % candidate region
      % Top axis point of the region:
      top = Points(Ind(i+j-1,1):Ind(i+j,2));
      Top = average(P(top,:));
      % Axis of the cylinder:
      Axis = Top-Bot;
      c0.axis = Axis/norm(Axis);
      % Compute the height along the axis:
      h = (P(RegC,:)-Bot)*c0.axis';
      minh = min(h);
      % Correct Bot to correspond to the real bottom
      if j == 0
        Bot = Bot+minh*c0.axis;
        c0.start = Bot;
        h = (P(RegC,:)-Bot)*c0.axis';
        minh = min(h);
      end
      if i+j >= nl
        ht = (Top-c0.start)*c0.axis';
        Top = Top+(max(h)-ht)*c0.axis;
      end
      % Compute the height of the Top:
      ht = (Top-c0.start)*c0.axis';
      Sec = h <= ht & h >= minh; % only points below the Top
      c0.length = ht-minh; % length of the region/cylinder
      % The region for the cylinder fitting:
      reg = RegC(Sec);
      Q0 = P(reg,:);

      %% Filter points and estimate radius
      if size(Q0,1) > 20
        [Keep,c0] = surface_coverage_filtering(Q0,c0,0.02,20);
        reg = reg(Keep);
        Q0 = Q0(Keep,:);
      else
        c0.radius = 0.01;
        c0.SurfCov = 0.05;
        c0.mad = 0.01;
        c0.conv = 1;
        c0.rel = 1;
      end

      %% Fit cylinder
      if size(Q0,1) > 9
        if i >= nl && t == 0
          c = least_squares_cylinder(Q0,c0);
        elseif i >= nl && t > 0
          h = (Q0-CylTop)*c0.axis';
          I = h >= 0;
          Q = Q0(I,:); % the section
          reg = reg(I);
          n2 = size(Q,1);     n1 = nnz(~I);
          if n2 > 9 && n1 > 5
            Q0 = [Q0(~I,:); Q]; % the point cloud for cylinder fitting
            W = [1/3*ones(n2,1); 2/3*ones(n1,1)]; % the weights
            c = least_squares_cylinder(Q0,c0,W,Q);
          else
            c = least_squares_cylinder(Q0,c0);
          end
        elseif t == 0
          top = Points(Ind(i+j-3,1):Ind(i+j-2,2));
          Top = average(P(top,:)); % Top axis point of the region
          ht = (Top-Bot)*c0.axis';
          h = (Q0-Bot)*c0.axis';
          I = h <= ht;
          Q = Q0(I,:); % the section
          reg = reg(I);
          n2 = size(Q,1);     n3 = nnz(~I);
          if n2 > 9 && n3 > 5
            Q0 = [Q; Q0(~I,:)]; % the point cloud for cylinder fitting
            W = [2/3*ones(n2,1); 1/3*ones(n3,1)]; % the weights
            c = least_squares_cylinder(Q0,c0,W,Q);
          else
            c = least_squares_cylinder(Q0,c0);
          end
        else
          top = Points(Ind(i+j-3,1):Ind(i+j-2,2));
          Top = average(P(top,:)); % Top axis point of the region
          ht = (Top-CylTop)*c0.axis';
          h = (Q0-CylTop)*c0.axis';
          I1 = h < 0; % the bottom
          I2 = h >= 0 & h <= ht; % the section
          I3 = h > ht; % the top
          Q = Q0(I2,:);
          reg = reg(I2);
          n1 = nnz(I1);   n2 = size(Q,1);     n3 = nnz(I3);
          if n2 > 9
            Q0 = [Q0(I1,:); Q; Q0(I3,:)];
            W = [1/4*ones(n1,1); 2/4*ones(n2,1); 1/4*ones(n3,1)];
            c = least_squares_cylinder(Q0,c0,W,Q);
          else
            c = c0;
            c.rel = 0;
          end
        end

        if c.conv == 0
          c = c0;
          c.rel = 0;
        end
        if c.SurfCov < 0.2
          c.rel = 0;
        end
      else
        c = c0;
        c.rel = 0;
      end

      % Collect fit data
      data(j+1,:) = [c.rel c.conv c.SurfCov c.length/c.radius];
      cyls{j+1} = c;
      regs{j+1} = reg;
      j = j+1;
      % If reasonable cylinder fitted, then stop fitting new ones
      % (but always fit at least three cylinders)
      RL = c.length/c.radius; % relative length of the cylinder
      if again && c.rel && c.conv && RL > 2
        if si == 1 && c.SurfCov > 0.7
          again = false;
        elseif si > 1 && c.SurfCov > 0.5
          again = false;
        end
      end
    end

    %% Select the best of the fitted cylinders
    % based on maximum surface coverage
    OKfit = data(1:j,1) & data(1:j,2) & data(1:j,4) > 1.5;

    J = (1:1:j)';
    t = t+1;
    if any(OKfit)
      J = J(OKfit);
    end
    [~,I] = max(data(J,3)-0.01*data(J,4));
    J = J(I);
    c = cyls{J};

    %% Update the indexes of the layers for the next region:
    CylTop = c.start+c.length*c.axis;
    i0 = i0+1;
    bot = Points(Ind(i0,1):Ind(i0+1,2));
    Bot = average(P(bot,:)); % Bottom axis point of the region
    h = (Bot-CylTop)*c.axis';
    i00 = i0;
    while i0+1 < nl && i0 < i00+5 && h < -c.radius/3
      i0 = i0+1;
      bot = Points(Ind(i0,1):Ind(i0+1,2));
      Bot = average(P(bot,:)); % Bottom axis point of the region
      h = (Bot-CylTop)*c.axis';
    end
    i = i0+5;
    i = min(i,nl);

    %% If the next section is very short part of the end of the branch
    % then simply increase the length of the current cylinder
    if nl-i0+2 < 4
      reg = Points(Ind(nl-5,1):Ind(nl,2));
      Q0 = P(reg,:);
      ht = (c.start+c.length*c.axis)*c.axis';
      h = Q0*c.axis';
      maxh = max(h);
      if maxh > ht
        c.length = c.length+(maxh-ht);
      end
      i0 = nl;
    end
    Reg{t} = regs{J};

    if t == 1
      cyl = c;
      names = fieldnames(cyl);
      n = max(size(names));
    else
      for k = 1:n
        cyl.(names{k}) = [cyl.(names{k}); c.(names{k})];
      end
    end

    %% compute cylinder top for the definition of the next section
    CylTop = c.start+c.length*c.axis;
  end
  Reg = Reg(1:t);

else
  %% Define a region for small segments
  Q0 = P(Points,:);
  if size(Q0,1) > 10
    %% Define the direction
    bot = Points(Ind(1,1):Ind(1,2));
    Bot = average(P(bot,:));
    top = Points(Ind(nl,1):Ind(nl,2));
    Top = average(P(top,:));
    Axis = Top-Bot;
    c0.axis = Axis/norm(Axis);
    h = Q0*c0.axis';
    c0.length = max(h)-min(h);
    hpoint = Bot*c0.axis';
    c0.start = Bot-(hpoint-min(h))*c0.axis;

    %% Define other outputs
    [Keep,c0] = surface_coverage_filtering(Q0,c0,0.02,20);
    Reg = cell(1,1);
    Reg{1} = Points(Keep);
    Q0 = Q0(Keep,:);
    cyl = least_squares_cylinder(Q0,c0);
    if ~cyl.conv || ~cyl.rel
      cyl = c0;
    end
    t = 1;
  else
    cyl = 0;
    t = 0;
  end
end
% Define Reg as coordinates
for i = 1:t
  Reg{i} = P(Reg{i},:);
end
Reg = Reg(1:t);
% End of function
end


function [PC,cyl,added] = parent_cylinder(SPar,SChi,CiS,cylinder,cyl,si)

% Finds the parent cylinder from the possible parent segment.
% Does this by checking if the axis of the cylinder, if continued, will
% cross the nearby cylinders in the parent segment.
% Adjust the cylinder so that it starts from the surface of its parent.

rad = cyl.radius;
len = cyl.length;
sta = cyl.start;
axe = cyl.axis;

% PC     Parent cylinder
nc = numel(rad);
added = false;
if SPar(si) > 0 % parent segment exists, find the parent cylinder
  s = SPar(si);
  PC = CiS{s}; % the cylinders in the parent segment
  % select the closest cylinders for closer examination
  if length(PC) > 1
    D = mat_vec_subtraction(-cylinder.start(PC,:),-sta(1,:));
    d = sum(D.*D,2);
    [~,I] = sort(d);
    if length(PC) > 3
      I = I(1:4);
    end
    pc = PC(I);
    ParentFound = false;
  elseif length(PC) == 1
    ParentFound = true;
  else
    PC = zeros(0,1);
    ParentFound = true;
  end

  %% Check possible crossing points
  if ~ParentFound
    pc0 = pc;
    n = length(pc);
    % Calculate the possible crossing points of the cylinder axis, when
    % extended, on the surfaces of the parent candidate cylinders
    x = zeros(n,2);  % how much the starting point has to move to cross
    h = zeros(n,2);  % the crossing point height in the parent
    Axe = cylinder.axis(pc,:);
    Sta = cylinder.start(pc,:);
    for j = 1:n
      % Crossing points solved from a quadratic equation
      A = axe(1,:)-(axe(1,:)*Axe(j,:)')*Axe(j,:);
      B = sta(1,:)-Sta(j,:)-(sta(1,:)*Axe(j,:)')*Axe(j,:)...
        +(Sta(j,:)*Axe(j,:)')*Axe(j,:);
      e = A*A';
      f = 2*A*B';
      g = B*B'-cylinder.radius(pc(j))^2;
      di = sqrt(f^2 - 4*e*g);  % the discriminant
      s1 = (-f + di)/(2*e);
      % how much the starting point must be moved to cross:
      s2 = (-f - di)/(2*e);
      if isreal(s1) %% cylinders can cross
        % the heights of the crossing points
        x(j,:) = [s1 s2];
        h(j,1) = sta(1,:)*Axe(j,:)'+x(j,1)*axe(1,:)*Axe(j,:)'-...
          Sta(j,:)*Axe(j,:)';
        h(j,2) = sta(1,:)*Axe(j,:)'+x(j,2)*axe(1,:)*Axe(j,:)'-...
          Sta(j,:)*Axe(j,:)';
      end
    end

    %% Extend to crossing point in the (extended) parent
    I = x(:,1) ~= 0; % Select only candidates with crossing points
    pc = pc0(I);    x = x(I,:);     h = h(I,:);
    j = 1;      n = nnz(I);
    X = zeros(n,3); %
    Len = cylinder.length(pc);
    while j <= n && ~ParentFound
      if x(j,1) > 0 && x(j,2) < 0
        % sp inside the parent and crosses its surface
        if h(j,1) >= 0 && h(j,1) <= Len(j) && len(1)-x(j,1) > 0
          PC = pc(j);
          sta(1,:) = sta(1,:)+x(j,1)*axe(1,:);
          len(1) = len(1)-x(j,1);
          ParentFound = true;
        elseif len(1)-x(j,1) > 0
          if h(j,1) < 0
            X(j,:) = [x(j,1) abs(h(j,1)) 0];
          else
            X(j,:) = [x(j,1) h(j,1)-Len(j) 0];
          end
        else
          X(j,:) = [x(j,1) h(j,1) 1];
        end
      elseif x(j,1) < 0 && x(j,2) > 0 && len(1)-x(j,2) > 0
        % sp inside the parent and crosses its surface
        if h(j,2) >= 0 && h(j,2) <= Len(j) && len(1)-x(j,2) > 0
          PC = pc(j);
          sta(1,:) = sta(1,:)+x(j,2)*axe(1,:);
          len(1) = len(1)-x(j,2);
          ParentFound = true;
        elseif len(1)-x(j,2) > 0
          if h(j,2) < 0
            X(j,:) = [x(j,2) abs(h(j,2)) 0];
          else
            X(j,:) = [x(j,2) h(j,2)-Len(j) 0];
          end
        else
          X(j,:) = [x(j,2) h(j,2) 1];
        end
      elseif x(j,1) < 0 && x(j,2) < 0 && x(j,2) < x(j,1) && len(1)-x(j,1) > 0
        % sp outside the parent and crosses its surface when extended
        % backwards
        if h(j,1) >= 0 && h(j,1) <= Len(j) && len(1)-x(j,1) > 0
          PC = pc(j);
          sta(1,:) = sta(1,:)+x(j,1)*axe(1,:);
          len(1) = len(1)-x(j,1);
          ParentFound = true;
        elseif len(1)-x(j,1) > 0
          if h(j,1) < 0
            X(j,:) = [x(j,1) abs(h(j,1)) 0];
          else
            X(j,:) = [x(j,1) h(j,1)-Len(j) 0];
          end
        else
          X(j,:) = [x(j,1) h(j,1) 1];
        end
      elseif x(j,1) < 0 && x(j,2) < 0 && x(j,2) > x(j,1) && len(1)-x(j,2) > 0
        % sp outside the parent and crosses its surface when extended
        % backwards
        if h(j,2) >= 0 && h(j,2) <= Len(j) && len(1)-x(j,2) > 0
          PC = pc(j);
          sta(1,:) = sta(1,:)+x(j,2)*axe(1,:);
          len(1) = len(1)-x(j,2);
          ParentFound = true;
        elseif len(1)-x(j,2) > 0
          if h(j,2) < 0
            X(j,:) = [x(j,2) abs(h(j,2)) 0];
          else
            X(j,:) = [x(j,2) h(j,2)-Len(j) 0];
          end
        else
          X(j,:) = [x(j,2) h(j,2) 1];
        end
      elseif x(j,1) > 0 && x(j,2) > 0 && x(j,2) < x(j,1) && len(1)-x(j,1) > 0
        % sp outside the parent but crosses its surface when extended forward
        if h(j,1) >= 0 && h(j,1) <= Len(j) && len(1)-x(j,1) > 0
          PC = pc(j);
          sta(1,:) = sta(1,:)+x(j,1)*axe(1,:);
          len(1) = len(1)-x(j,1);
          ParentFound = true;
        elseif len(1)-x(j,1) > 0
          if h(j,1) < 0
            X(j,:) = [x(j,1) abs(h(j,1)) 0];
          else
            X(j,:) = [x(j,1) h(j,1)-Len(j) 0];
          end
        else
          X(j,:) = [x(j,1) h(j,1) 1];
        end
      elseif x(j,1) > 0 && x(j,2) > 0 && x(j,2) > x(j,1) && len(1)-x(j,2) > 0
        % sp outside the parent and crosses its surface when extended forward
        if h(j,2) >= 0 && h(j,2) <= Len(j) && len(1)-x(j,2) > 0
          PC = pc(j);
          sta(1,:) = sta(1,:)+x(j,2)*axe(1,:);
          len(1) = len(1)-x(j,2);
          ParentFound = true;
        elseif len(1)-x(j,2) > 0
          if h(j,1) < 0
            X(j,:) = [x(j,2) abs(h(j,2)) 0];
          else
            X(j,:) = [x(j,2) h(j,2)-Len(j) 0];
          end
        else
          X(j,:) = [x(j,2) h(j,2) 1];
        end
      end
      j = j+1;
    end

    if ~ParentFound && n > 0
      [H,I] = min(X(:,2));
      X = X(I,:);
      if X(3) == 0 && H < 0.1*Len(I)
        PC = pc(I);
        sta(1,:) = sta(1,:)+X(1)*axe(1,:);
        len(1) = len(1)-X(1);
        ParentFound = true;
      else
        PC = pc(I);

        if nc > 1 && X(1) <= rad(1) && abs(X(2)) <= 1.25*cylinder.length(PC)
          % Remove the first cylinder and adjust the second
          S = sta(1,:)+X(1)*axe(1,:);
          V = sta(2,:)+len(2)*axe(2,:)-S;
          len(2) = norm(V);         len = len(2:nc);
          axe(2,:) = V/norm(V);      axe = axe(2:nc,:);
          sta(2,:) = S;            sta = sta(2:nc,:);
          rad = rad(2:nc);
          cyl.mad = cyl.mad(2:nc);
          cyl.SurfCov = cyl.SurfCov(2:nc);
          nc = nc-1;
          ParentFound = true;
        elseif nc > 1
          % Remove the first cylinder
          sta = sta(2:nc,:);    len = len(2:nc);
          axe = axe(2:nc,:);        rad = rad(2:nc);
          cyl.mad = cyl.mad(2:nc);
          cyl.SurfCov = cyl.SurfCov(2:nc);
          nc = nc-1;
        elseif isempty(SChi{si})
          % Remove the cylinder
          nc = 0;
          PC = zeros(0,1);
          ParentFound = true;
          rad = zeros(0,1);
        elseif X(1) <= rad(1) && abs(X(2)) <= 1.5*cylinder.length(PC)
          % Adjust the cylinder
          sta(1,:) = sta(1,:)+X(1)*axe(1,:);
          len(1) = abs(X(1));
          ParentFound = true;
        end
      end
    end

    if ~ParentFound
      % The parent is the cylinder in the parent segment whose axis
      % line is the closest to the axis line of the first cylinder
      % Or the parent cylinder is the one whose base, when connected
      % to the first cylinder is the most parallel.
      % Add new cylinder
      pc = pc0;

      [Dist,~,DistOnLines] = distances_between_lines(...
        sta(1,:),axe(1,:),cylinder.start(pc,:),cylinder.axis(pc,:));

      I = DistOnLines >= 0;
      J = DistOnLines <= cylinder.length(pc);
      I = I&J;
      if ~any(I)
        I = DistOnLines >= -0.2*cylinder.length(pc);
        J = DistOnLines <= 1.2*cylinder.length(pc);
        I = I&J;
      end
      if any(I)
        pc = pc(I);     Dist = Dist(I);     DistOnLines = DistOnLines(I);
        [~,I] = min(Dist);
        DistOnLines = DistOnLines(I);       PC = pc(I);
        Q = cylinder.start(PC,:)+DistOnLines*cylinder.axis(PC,:);
        V = sta(1,:)-Q;      L = norm(V);        V = V/L;
        a = acos(V*cylinder.axis(PC,:)');
        h = sin(a)*L;
        S = Q+cylinder.radius(PC)/h*L*V;
        L = (h-cylinder.radius(PC))/h*L;
        if L > 0.01 && L/len(1) > 0.2
          nc = nc+1;
          sta = [S; sta];   rad = [rad(1); rad];
          axe = [V; axe];       len = [L; len];
          cyl.mad = [cyl.mad(1); cyl.mad];
          cyl.SurfCov = [cyl.SurfCov(1); cyl.SurfCov];
          cyl.rel = [cyl.rel(1); cyl.rel];
          cyl.conv = [cyl.conv(1); cyl.conv];
          added = true;
        end
      else
        V = -mat_vec_subtraction(cylinder.start(pc,:),sta(1,:));
        L0 = sqrt(sum(V.*V,2));
        V = [V(:,1)./L0 V(:,2)./L0 V(:,3)./L0];
        A = V*axe(1,:)';
        [A,I] = max(A);
        L1 = L0(I);       PC = pc(I);     V = V(I,:);
        a = acos(V*cylinder.axis(PC,:)');
        h = sin(a)*L1;
        S = cylinder.start(PC,:)+cylinder.radius(PC)/h*L1*V;
        L = (h-cylinder.radius(PC))/h*L1;
        if L > 0.01 && L/len(1) > 0.2
          nc = nc+1;
          sta = [S; sta];   rad = [rad(1); rad];
          axe = [V; axe];   len = [L; len];
          cyl.mad = [cyl.mad(1); cyl.mad];
          cyl.SurfCov = [cyl.SurfCov(1); cyl.SurfCov];
          cyl.rel = [cyl.rel(1); cyl.rel];
          cyl.conv = [cyl.conv(1); cyl.conv];
          added = true;
        end
      end
    end
  end
else
  % no parent segment exists
  PC = zeros(0,1);
end

% define the output
cyl.radius = rad(1:nc);     cyl.length = len(1:nc,:);
cyl.start = sta(1:nc,:);    cyl.axis = axe(1:nc,:);
cyl.mad = cyl.mad(1:nc);    cyl.SurfCov = cyl.SurfCov(1:nc);
cyl.conv = cyl.conv(1:nc);  cyl.rel = cyl.rel(1:nc);
% End of function
end


function cyl = adjustments(cyl,parcyl,inputs,Regs)

nc = size(cyl.radius,1);
Mod = false(nc,1); % cylinders modified
SC = cyl.SurfCov;

%% Determine the maximum and the minimum radius
% The maximum based on parent branch
if ~isempty(parcyl.radius)
  MaxR = 0.95*parcyl.radius;
  MaxR = max(MaxR,inputs.MinCylRad);
else
  % use the maximum from the bottom cylinders
  a = min(3,nc);
  MaxR = 1.25*max(cyl.radius(1:a));
end
MinR = min(cyl.radius(SC > 0.7));
if ~isempty(MinR) && min(cyl.radius) < MinR/2
  MinR = min(cyl.radius(SC > 0.4));
elseif isempty(MinR)
  MinR = min(cyl.radius(SC > 0.4));
  if isempty(MinR)
    MinR = inputs.MinCylRad;
  end
end

%% Check maximum and minimum radii
I = cyl.radius < MinR;
cyl.radius(I) = MinR;
Mod(I) = true;
if inputs.ParentCor || nc <= 3
  I = (cyl.radius > MaxR & SC < 0.7) | (cyl.radius > 1.2*MaxR);
  cyl.radius(I) = MaxR;
  Mod(I) = true;
  % For short branches modify with more restrictions
  if nc <= 3
    I = (cyl.radius > 0.75*MaxR & SC < 0.7);
    if any(I)
      r = max(SC(I)/0.7.*cyl.radius(I),MinR);
      cyl.radius(I) = r;
      Mod(I) = true;
    end
  end
end

%% Use taper correction to modify radius of too small and large cylinders
% Adjust radii if a small SurfCov and high SurfCov in the previous and
% following cylinders
for i = 2:nc-1
  if SC(i) < 0.7 && SC(i-1) >= 0.7 && SC(i+1) >= 0.7
    cyl.radius(i) = 0.5*(cyl.radius(i-1)+cyl.radius(i+1));
    Mod(i) = true;
  end
end

%% Use taper correction to modify radius of too small and large cylinders
if inputs.TaperCor
  if max(cyl.radius) < 0.001

    %% Adjust radii of thin branches to be linearly decreasing
    if nc > 2
      r = sort(cyl.radius);
      r = r(2:end-1);
      a = 2*mean(r);
      if a > max(r)
        a = min(0.01,max(r));
      end
      b = min(0.5*min(cyl.radius),0.001);
      cyl.radius = linspace(a,b,nc)';
    elseif nc > 1
      r = max(cyl.radius);
      cyl.radius = [r; 0.5*r];
    end
    Mod = true(nc,1);

  elseif nc > 4
    %% Parabola adjustment of maximum and minimum
    % Define parabola taper shape as maximum (and minimum) radii for
    % the cylinders with low surface coverage
    branchlen = sum(cyl.length(1:nc)); % branch length
    L = cyl.length/2+[0; cumsum(cyl.length(1:nc-1))];
    Taper = [L; branchlen];
    Taper(:,2) = [1.05*cyl.radius; MinR];
    sc = [SC; 1];

    % Least square fitting of parabola to "Taper":
    A = [sum(sc.*Taper(:,1).^4) sum(sc.*Taper(:,1).^2); ...
      sum(sc.*Taper(:,1).^2) sum(sc)];
    y = [sum(sc.*Taper(:,2).*Taper(:,1).^2); sum(sc.*Taper(:,2))];
    warning off
    x = A\y;
    warning on
    x(1) = min(x(1),-0.0001); % tapering from the base to the tip
    Ru = x(1)*L.^2+x(2); % upper bound parabola
    Ru( Ru < MinR ) = MinR;
    if max(Ru) > MaxR
      a = max(Ru);
      Ru = MaxR/a*Ru;
    end
    Rl = 0.75*Ru; % lower bound parabola
    Rl( Rl < MinR ) = MinR;

    % Modify radii based on parabola:
    % change values larger than the parabola-values when SC < 70%:
    I = cyl.radius > Ru & SC < 0.7;
    cyl.radius(I) = Ru(I)+(cyl.radius(I)-Ru(I)).*SC(I)/0.7;
    Mod(I) = true;
    % change values larger than the parabola-values when SC > 70% and
    % radius is over 33% larger than the parabola-value:
    I = cyl.radius > 1.333*Ru & SC >= 0.7;
    cyl.radius(I) = Ru(I)+(cyl.radius(I)-Ru(I)).*SC(I);
    Mod(I) = true;
    % change values smaller than the downscaled parabola-values:
    I = (cyl.radius < Rl & SC < 0.7) | (cyl.radius < 0.5*Rl);
    cyl.radius(I) = Rl(I);
    Mod(I) = true;

  else
    %% Adjust radii of short branches to be linearly decreasing
    R = cyl.radius;
    if nnz(SC >= 0.7) > 1
      a = max(R(SC >= 0.7));
      b = min(R(SC >= 0.7));
    elseif nnz(SC >= 0.7) == 1
      a = max(R(SC >= 0.7));
      b = min(R);
    else
      a = sum(R.*SC/sum(SC));
      b = min(R);
    end
    Ru = linspace(a,b,nc)';
    I = SC < 0.7 & ~Mod;
    cyl.radius(I) = Ru(I)+(R(I)-Ru(I)).*SC(I)/0.7;
    Mod(I) = true;

  end
end

%% Modify starting points by optimising them for given radius and axis
nr = size(Regs,1);
for i = 1:nc
  if Mod(i)
    if nr == nc
      Reg = Regs{i};
    elseif i > 1
      Reg = Regs{i-1};
    end
    if abs(cyl.radius(i)-cyl.radius0(i)) > 0.005 && ...
        (nr == nc || (nr < nc && i > 1))
      P = Reg-cyl.start(i,:);
      [U,V] = orthonormal_vectors(cyl.axis(i,:));
      P = P*[U V];
      cir = least_squares_circle_centre(P,[0 0],cyl.radius(i));
      if cir.conv && cir.rel
        cyl.start(i,:) = cyl.start(i,:)+cir.point(1)*U'+cir.point(2)*V';
        cyl.mad(i,1) = cir.mad;
        [~,V,h] = distances_to_line(Reg,cyl.axis(i,:),cyl.start(i,:));
        if min(h) < -0.001
          cyl.length(i) = max(h)-min(h);
          cyl.start(i,:) = cyl.start(i,:)+min(h)*cyl.axis(i,:);
          [~,V,h] = distances_to_line(Reg,cyl.axis(i,:),cyl.start(i,:));
        end
        a = max(0.02,0.2*cyl.radius(i));
        nl = ceil(cyl.length(i)/a);
        nl = max(nl,4);
        ns = ceil(2*pi*cyl.radius(i)/a);
        ns = max(ns,10);
        ns = min(ns,36);
        cyl.SurfCov(i,1) = surface_coverage2(...
          cyl.axis(i,:),cyl.length(i),V,h,nl,ns);
      end
    end
  end
end

%% Continuous branches
% Make cylinders properly "continuous" by moving the starting points
% Move the starting point to the plane defined by parent cylinder's top
if nc > 1
  for j = 2:nc
    U = cyl.start(j,:)-cyl.start(j-1,:)-cyl.length(j-1)*cyl.axis(j-1,:);
    if (norm(U) > 0.0001)
      % First define vector V and W which are orthogonal to the
      % cylinder axis N
      N = cyl.axis(j,:)';
      if norm(N) > 0
        [V,W] = orthonormal_vectors(N);
        % Now define the new starting point
        x = [N V W]\U';
        cyl.start(j,:) = cyl.start(j,:)-x(1)*N';
        if x(1) > 0
          cyl.length(j) = cyl.length(j)+x(1);
        elseif cyl.length(j)+x(1) > 0
          cyl.length(j) = cyl.length(j)+x(1);
        end
      end
    end
  end
end

%% Connect far away first cylinder to the parent
if ~isempty(parcyl.radius)
  [d,V,h,B] = distances_to_line(cyl.start(1,:),parcyl.axis,parcyl.start);
  d = d-parcyl.radius;
  if d > 0.001
    taper = cyl.start(1,:);
    E = taper+cyl.length(1)*cyl.axis(1,:);
    V = parcyl.radius*V/norm(V);
    if h >= 0 && h <= parcyl.length
      cyl.start(1,:) = parcyl.start+B+V;
    elseif h < 0
      cyl.start(1,:) = parcyl.start+V;
    else
      cyl.start(1,:) = parcyl.start+parcyl.length*parcyl.axis+V;
    end
    cyl.axis(1,:) = E-cyl.start(1,:);
    cyl.length(1) = norm(cyl.axis(1,:));
    cyl.axis(1,:) = cyl.axis(1,:)/cyl.length(1);
  end
end

% End of function
end
