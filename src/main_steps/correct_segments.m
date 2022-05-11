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

function segment = correct_segments(P,cover,segment,inputs,RemSmall,ModBases,AddChild)

% ---------------------------------------------------------------------
% CORRECT_SEGMENTS.M        Corrects the given segmentation.
%
% Version 2.0.2
% Latest update     2 May 2022
%
% Copyright (C) 2013-2022 Pasi Raumonen
% ---------------------------------------------------------------------

% First segments are modified by making them as long as possible. Here the
% stem and 1-st order branches are handled differently as there is also
% restriction to how "curved" they can be in the sense of ratio
% total_length/base_tip_distance. Then, optionally, small segments that
% are close to their parent and have no children are removed as unclear
% (are they part of the parent or real segments?).
% Then, optionally, the bases of branches are modified by
% expanding them into parent segment in order to remove ledges from the
% parent from locations of the branches.

% Inputs:
% P             Point cloud
% cover         Cover sets
% segment       Segments
% inputs        The input structure
% RemSmall      If True, small unclear segments are removed
% ModBase       If True, bases of the segments are modified
% AddChild      If True, the expanded (modified) base is added to the child segment.
%               If AddChild = false and ModBase = true, then the expanded part is
%               removed from both the child and the parent.
% Outputs:
% segment       Segments
% ---------------------------------------------------------------------

% Changes from version 2.0.1 to 2.0.2, 2 May 2022:
% 1) Added "if ~isempty(SegChildren)... " statement to the
%    "modify_topology" subfunction where next branch is selected based on 
%    the increasing branching order to prevent a rare bug 

% Changes from version 2.0.0 to 2.0.1, 2 Oct 2019:
% 1) Main function: added "if SPar(i,1) > 1"-statement to ModBase -->
%    NotAddChild

if nargin == 4
  RemSmall = true;
  ModBases = false;
elseif nargin == 5
  ModBases = false;
elseif nargin == 6
  AddChild = false;
end

Bal = cover.ball;
Segs = segment.segments;
SPar = segment.ParentSegment;
SChi = segment.ChildSegment;
Ce = P(cover.center,:);

%% Make stem and branches as long as possible
if RemSmall
  [Segs,SPar,SChi] = modify_topology(P,Ce,Bal,Segs,SPar,SChi,inputs.PatchDiam2Max);
else
  [Segs,SPar,SChi] = modify_topology(P,Ce,Bal,Segs,SPar,SChi,inputs.PatchDiam1);
end

%% Remove small child segments
if RemSmall
  [Segs,SPar,SChi] = remove_small(Ce,Segs,SPar,SChi);
end

% Check the consistency of empty vector sizes
ns = size(Segs,1);
for i = 1:ns
  if isempty(SChi{i})
    SChi{i} = zeros(0,1,'uint32');
  end
end

if ModBases
  %% Modify the base of the segments
  ns = size(Segs,1);
  base = cell(200,1);
  if AddChild
    % Add the expanded base to the child and remove it from the parent
    for i = 2:ns
      SegC = Segs{i};
      SegP = Segs{SPar(i,1)};
      [SegP,Base] = modify_parent(P,Bal,Ce,SegP,SegC,SPar(i,2),inputs.PatchDiam1,base);
      Segs{SPar(i,1)} = SegP;
      SegC{1} = Base;
      Segs{i} = SegC;
    end
  else
    % Only remove the expanded base from the parent
    for i = 2:ns
      if SPar(i,1) > 1
        SegC = Segs{i};
        SegP = Segs{SPar(i,1)};
        SegP = modify_parent(P,Bal,Ce,SegP,SegC,SPar(i,2),inputs.PatchDiam2Max,base);
        Segs{SPar(i,1)} = SegP;
      end
    end
  end
end
SPar = SPar(:,1);

% Modify the size and type of SChi and Segs, if necessary
ns = size(Segs,1);
for i = 1:ns
  C = SChi{i};
  if size(C,2) > size(C,1) && size(C,1) > 0
    SChi{i} = uint32(C');
  elseif size(C,1) == 0 || size(C,2) == 0
    SChi{i} = zeros(0,1,'uint32');
  else
    SChi{i} = uint32(C);
  end
  S = Segs{i};
  for j = 1:size(S,1)
    S{j} = uint32(S{j});
  end
  Segs{i} = S;
end
segment.segments = Segs;
segment.ParentSegment = SPar;
segment.ChildSegment = SChi;

%% Generate segment data for the points
np = size(P,1);
ns = size(Segs,1);
% Define for each point its segment
if ns <= 2^16
  SegmentOfPoint = zeros(np,1,'uint16');
else
  SegmentOfPoint = zeros(np,1,'uint32');
end
for i = 1:ns
  S = Segs{i};
  S = vertcat(S{:});
  SegmentOfPoint(vertcat(Bal{S})) = i;
end
segment.SegmentOfPoint = SegmentOfPoint;
% Define the indexes of the segments up to 3rd-order
C = SChi{1};
segment.branch1indexes = C;
if ~isempty(C)
  C = vertcat(SChi{C});
  segment.branch2indexes = C;
  if ~isempty(C)
    C = vertcat(SChi{C});
    segment.branch3indexes = C;
  else
    segment.branch3indexes = zeros(0,1);
  end
else
  segment.branch2indexes = zeros(0,1);
  segment.branch3indexes = zeros(0,1);
end

end % End of main function


function StemTop = search_stem_top(P,Ce,Bal,Segs,SPar,dmin)

% Search the stem's top segment such that the resulting stem
% 1) is one the highest segments (goes to the top of the tree)
% 2) is horizontally close to the bottom of the stem (goes straigth up)
% 3) has length close to the distance between its bottom and top (is not too curved)
nseg = size(Segs,1);
SegHeight = zeros(nseg,1); % heights of the tips of the segments
HorDist = zeros(nseg,1); % horizontal distances of the tips from stem's center
s = Segs{1}{1};
StemCen = average(Ce(s,:)); % center (x,y) of stem base
for i = 1:nseg
  S = Segs{i}{end}(1);
  SegHeight(i) = Ce(S,3);
  HorDist(i) = norm(Ce(S,1:2)-StemCen(1:2));
end
Top = max(SegHeight); % the height of the highest tip
HeiDist = Top-SegHeight; % the height difference to "Top"
Dist = sqrt((HorDist.^2+HeiDist.^2)); % Distance to the top
LenDisRatio = 2;
SearchDist = 0.5;
MaxLenDisRatio = 1.05; % the maximum acceptable length/distance ratio of segments
SubSegs = zeros(100,1); % Segments to be combined to form the stem
while LenDisRatio > MaxLenDisRatio
  StemTops = (1:1:nseg)';
  I = Dist < SearchDist; % only segments with distance to the top < 0.5m
  while ~any(I)
    SearchDist = SearchDist+0.5;
    I = Dist < SearchDist;
  end
  StemTops = StemTops(I);

  % Define i-1 alternative stems from StemTops
  n = length(StemTops);
  Stems = cell(n,1);
  Segment = cell(3000,1);
  for j = 1:n
    Seg = Segs{1};
    spar = SPar;
    if StemTops(j) ~= 1
      % Tip point was not in the current segment, modify segments
      SubSegs(1) = StemTops(j);
      nsegs = 1;
      segment = StemTops(j);
      while segment ~= 1
        segment = SPar(segment,1);
        nsegs = nsegs+1;
        SubSegs(nsegs) = segment;
      end
      % Modify stem
      a = size(Seg,1);
      Segment(1:a) = Seg;
      a = a+1;
      for i = 1:nsegs-2
        I = SubSegs(nsegs-i); % segment to be combined to the first segment
        J = SubSegs(nsegs-i-1); % above segment's child to be combined next
        SP = spar(I,2);  % layer index of the child in the parent
        SegC = Segs{I};
        sp = spar(J,2);  % layer index of the child's child in the child
        if SP >= a-2 % Use the whole parent
          Segment(a:a+sp-1) = SegC(1:sp);
          spar(J,2) = a+sp-1;
          a = a+sp;
        else % Use only bottom part of the parent
          Segment(SP+1:SP+sp) = SegC(1:sp);
          a = SP+sp+1;
          spar(J,2) = SP+sp;
        end
        SubSegs(nsegs-i) = 1;
      end

      % Combine the last segment to the branch
      I = SubSegs(1);
      SP = spar(I,2);
      SegC = Segs{I};
      nc = size(SegC,1);
      if SP >= a-2 % Use the whole parent
        Segment(a:a+nc-1) = SegC;
        a = a+nc-1;
      else % divide the parent segment into two parts
        Segment(SP+1:SP+nc) = SegC;
        a = SP+nc;
      end
      Stems{j,1} = Segment(1:a);
    else
      Stems{j,1} = Seg;
    end

  end

  % Calculate the lengths of the candidate stems
  N = ceil(0.5/dmin/1.4); % number of layers used for linear length approximation
  Lengths = zeros(n,1);
  Heights = zeros(n,1);
  for i = 1:n
    Seg = Stems{i,1};
    ns = size(Seg,1);
    if ceil(ns/N) > floor(ns/N)
      m = ceil(ns/N);
    else
      m = ceil(ns/N)+1;
    end
    Nodes = zeros(m,3);
    for j = 1:m
      I = (j-1)*N+1;
      if I > ns
        I = ns;
      end
      S = Seg{I};
      if length(S) > 1
        Nodes(j,:) = average(Ce(S,:));
      else
        S = Bal{S};
        Nodes(j,:) = average(P(S,:));
      end
    end
    V = Nodes(2:end,:)-Nodes(1:end-1,:);
    Lengths(i) = sum(sqrt(sum(V.*V,2)));
    V = Nodes(end,:)-Nodes(1,:);
    Heights(i) = norm(V);
  end

  LenDisRatio = Lengths./Heights;
  [LenDisRatio,I] = min(LenDisRatio);
  StemTop = StemTops(I);
  SearchDist = SearchDist+1;
  if SearchDist > 3
    MaxLenDisRatio = 1.1;
    if SearchDist > 5
      MaxLenDisRatio = 1.15;
      if SearchDist > 7
        MaxLenDisRatio = 5;
      end
    end
  end
end

end % End subfunction


function BranchTop = search_branch_top(P,Ce,Bal,Segs,SPar,SChi,dmin,BI)

% Search the end segment for branch such that the resulting branch
% 1) has length close to the distance between its bottom and top
% 2) has distance close to the farthest segment end

% Inputs
% BI    Branch (segment) index

% Outputs
% BranchTop     The index of the segment forming the tip of the branch
%                   originating from the base of the given segment BI

% Define all the sub-segments of the given segments
ns = size(Segs,1);
Segments = zeros(ns,1); % the given segment and its sub-segments
Segments(1) = BI;
t = 2;
C = SChi{BI};
while ~isempty(C)
  n = length(C);
  Segments(t:t+n-1) = C;
  C = vertcat(SChi{C});
  t = t+n;
end
if t > 2
  t = t-n;
end
Segments = Segments(1:t);

% Determine linear distances from the segment tips to the base of the given
% segment
LinearDist = zeros(t,1); % linear distances from the
Seg = Segs{Segments(1)};
BranchBase = average(Ce(Seg{1},:)); % center of branch's base
for i = 1:t
  Seg = Segs{Segments(i)};
  C = average(Ce(Seg{end},:)); % tip
  LinearDist(i) = norm(C-BranchBase);
end
LinearDist = LinearDist(1:t);

% Sort the segments according their linear distance, from longest to
% shortest
[LinearDist,I] = sort(LinearDist,'descend');
Segments = Segments(I);

% Define alternative branches from Segments
Branches = cell(t,1); % the alternative segments as cell layers
SubSegs = zeros(100,1); % Segments to be combined
Segment = cell(3000,1);
for j = 1:t
  Seg = Segs{BI};
  spar = SPar;
  if Segments(j) ~= BI
    % Tip point was not in the current segment, modify segments
    SubSegs(1) = Segments(j);
    k = 1;
    S = Segments(j);
    while S ~= BI
      S = SPar(S,1);
      k = k+1;
      SubSegs(k) = S;
    end
    % Modify branch
    a = size(Seg,1);
    Segment(1:a) = Seg;
    a = a+1;
    for i = 1:k-2
      I = SubSegs(k-i); % segment to be combined to the first segment
      J = SubSegs(k-i-1); % above segment's child to be combined next
      SP = spar(I,2);  % layer index of the child in the parent
      SegC = Segs{I};
      sp = spar(J,2);  % layer index of the child's child in the child
      if SP >= a-2 % Use the whole parent
        Segment(a:a+sp-1) = SegC(1:sp);
        spar(J,2) = a+sp-1;
        a = a+sp;
      else % Use only bottom part of the parent
        Segment(SP+1:SP+sp) = SegC(1:sp);
        a = SP+sp+1;
        spar(J,2) = SP+sp;
      end
      SubSegs(k-i) = 1;
    end

    % Combine the last segment to the branch
    I = SubSegs(1);
    SP = spar(I,2);
    SegC = Segs{I};
    L = size(SegC,1);
    if SP >= a-2 % Use the whole parent
      Segment(a:a+L-1) = SegC;
      a = a+L-1;
    else % divide the parent segment into two parts
      Segment(SP+1:SP+L) = SegC;
      a = SP+L;
    end
    Branches{j,1} = Segment(1:a);
  else
    Branches{j,1} = Seg;
  end

end

% Calculate the lengths of the candidate branches. Stop, if possible, when
% the ratio length/linear distance is less 1.2 (branch is quite straight)
N = ceil(0.25/dmin/1.4); % number of layers used for linear length approximation
i = 1; % running index for while loop
Continue = true; % continue while loop as long as "Continue" is true
Lengths = zeros(t,1);  % linear lengths of the branches
while i <= t && Continue
  % Approximate the length with line segments connecting nodes along
  % the segment
  Seg = Branches{i,1};
  ns = size(Seg,1);
  if ceil(ns/N) > floor(ns/N)
    m = ceil(ns/N);
  else
    m = ceil(ns/N)+1;
  end
  Nodes = zeros(m,3);
  for j = 1:m
    I = (j-1)*N+1;
    if I > ns
      I = ns;
    end
    S = Seg{I};
    if length(S) > 1
      Nodes(j,:) = average(Ce(S,:));
    else
      S = Bal{S};
      Nodes(j,:) = average(P(S,:));
    end
  end
  V = Nodes(2:end,:)-Nodes(1:end-1,:); % line segments
  Lengths(i) = sum(sqrt(sum(V.*V,2)));

  % Continue as long as the length is less than 20% longer than the linear dist.
  % and the linear distance is over 75% of the maximum
  if Lengths(i)/LinearDist(i) < 1.20 && LinearDist(i) > 0.75*LinearDist(1)
    Continue = false;
    BranchTop = Segments(i);
  end
  i = i+1;
end

% If no suitable segment was found, try first with less strict conditions,
% and if that does not work, then select the one with the largest linear distance
if Continue
  L = Lengths./LinearDist;
  i = 1;
  while i <= t && L(i) > 1.4 && LinearDist(i) > 0.75*LinearDist(1)
    i = i+1;
  end
  if i <= t
    BranchTop = Segments(i);
  else
    BranchTop = Segments(1);
  end
end

end % End subfunction


function [Segs,SPar,SChi] = modify_topology(P,Ce,Bal,Segs,SPar,SChi,dmin)

% Make stem and branches as long as possible
ns = size(Segs,1);
Fal = false(2*ns,1);
nc = ceil(ns/5);
SubSegments = zeros(nc,1); % for searching sub-segments
SegInd = 1; % the segment under modification
UnMod = true(ns,1);
UnMod(SegInd) = false;
BranchOrder = 0;
ChildSegInd = 1; % index of the child segments under modification
while any(UnMod)
  ChildSegs = SChi{SegInd}; % child segments of the segment under modification
  if size(ChildSegs,1) < size(ChildSegs,2)
    ChildSegs = ChildSegs';
    SChi{SegInd} = ChildSegs;
  end

  if ~isempty(Segs(SegInd)) && ~isempty(ChildSegs)

    if SegInd > 1 && BranchOrder > 1 % 2nd-order and higher branches
      % Search the tip of the sub-branches with biggest linear
      % distance from the current branch's base
      SubSegments(1) = SegInd;
      NSubSegs = 2;
      while ~isempty(ChildSegs)
        n = length(ChildSegs);
        SubSegments(NSubSegs:NSubSegs+n-1) = ChildSegs;
        ChildSegs = vertcat(SChi{ChildSegs});
        NSubSegs = NSubSegs+n;
      end
      if NSubSegs > 2
        NSubSegs = NSubSegs-n;
      end

      % Find tip-points
      Top = zeros(NSubSegs,3);
      for i = 1:NSubSegs
        Top(i,:) = Ce(Segs{SubSegments(i)}{end}(1),:);
      end

      % Define bottom of the branch
      BotLayer = Segs{SegInd}{1};
      Bottom = average(Ce(BotLayer,:));

      % End segment is the segment whose tip has greatest distance to
      % the bottom of the branch
      V = mat_vec_subtraction(Top,Bottom);
      d = sum(V.*V,2);
      [~,I] = max(d);
      TipSeg = SubSegments(I(1));

    elseif SegInd > 1 && BranchOrder <= 1 % first order branches

      TipSeg = search_branch_top(P,Ce,Bal,Segs,SPar,SChi,dmin,SegInd);

    else % Stem

      TipSeg = search_stem_top(P,Ce,Bal,Segs,SPar,dmin);

    end

    if TipSeg ~= SegInd
      % Tip point was not in the current segment, modify segments
      SubSegments(1) = TipSeg;
      NSubSegs = 1;
      while TipSeg ~= SegInd
        TipSeg = SPar(TipSeg,1);
        NSubSegs = NSubSegs+1;
        SubSegments(NSubSegs) = TipSeg;
      end

      % refine branch
      for i = 1:NSubSegs-2
        I = SubSegments(NSubSegs-i); % segment to be combined to the first segment
        J = SubSegments(NSubSegs-i-1); % above segment's child to be combined next
        SP = SPar(I,2);  % layer index of the child in the parent
        SegP = Segs{SegInd};
        SegC = Segs{I};
        N = size(SegP,1);
        sp = SPar(J,2);  % layer index of the child's child in the child
        if SP >= N-2 % Use the whole parent
          Segs{SegInd} = [SegP; SegC(1:sp)];
          if sp < size(SegC,1) % use only part of the child segment
            Segs{I} = SegC(sp+1:end);
            SPar(I,2) = N+sp;

            ChildSegs = SChi{I};
            K = SPar(ChildSegs,2) <= sp;
            c = ChildSegs(~K);
            SChi{I} = c;
            SPar(c,2) = SPar(c,2)-sp;
            ChildSegs = ChildSegs(K);
            SChi{SegInd} = [SChi{SegInd}; ChildSegs];
            SPar(ChildSegs,1) = SegInd;
            SPar(ChildSegs,2) = N+SPar(ChildSegs,2);

          else % use the whole child segment
            Segs{I} = cell(0,1);
            SPar(I,1) = 0;
            UnMod(I) = false;

            ChildSegs = SChi{I};
            SChi{I} = zeros(0,1);
            c = set_difference(SChi{SegInd},I,Fal);
            SChi{SegInd} = [c; ChildSegs];
            SPar(ChildSegs,1) = SegInd;
            SPar(ChildSegs,2) = N+SPar(ChildSegs,2);

          end

          SubSegments(NSubSegs-i) = SegInd;
        else % divide the parent segment into two parts
          ns = ns+1;
          Segs{ns} = SegP(SP+1:end); % the top part of the parent forms a new segment
          SPar(ns,1) = SegInd;
          SPar(ns,2) = SP;
          UnMod(ns) = true;

          Segs{SegInd} = [SegP(1:SP); SegC(1:sp)];

          ChildSegs = SChi{SegInd};
          if size(ChildSegs,1) < size(ChildSegs,2)
            ChildSegs = ChildSegs';
          end
          K = SPar(ChildSegs,2) > SP;
          SChi{SegInd} = ChildSegs(~K);
          ChildSegs = ChildSegs(K);
          SChi{ns} = ChildSegs;
          SPar(ChildSegs,1) = ns;
          SPar(ChildSegs,2) = SPar(ChildSegs,2)-SP;
          SChi{SegInd} = [SChi{SegInd}; ns];
          if sp < size(SegC,1) % use only part of the child segment
            Segs{I} = SegC(sp+1:end);
            SPar(I,2) = SP+sp;

            ChildSegs = SChi{I};
            K = SPar(ChildSegs,2) <= sp;
            SChi{I} = ChildSegs(~K);
            SPar(ChildSegs(~K),2) = SPar(ChildSegs(~K),2)-sp;
            ChildSegs = ChildSegs(K);
            SChi{SegInd} = [SChi{SegInd}; ChildSegs];
            SPar(ChildSegs,1) = SegInd;
            SPar(ChildSegs,2) = SP+SPar(ChildSegs,2);

          else % use the whole child segment
            Segs{I} = cell(0,1);
            SPar(I,1) = 0;
            UnMod(I) = false;

            ChildSegs = SChi{I};
            c = set_difference(SChi{SegInd},I,Fal);
            SChi{SegInd} = [c; ChildSegs];
            SPar(ChildSegs,1) = SegInd;
            SPar(ChildSegs,2) = SP+SPar(ChildSegs,2);

          end
          SubSegments(NSubSegs-i) = SegInd;
        end

      end

      % Combine the last segment to the branch
      I = SubSegments(1);
      SP = SPar(I,2);
      SegP = Segs{SegInd};
      SegC = Segs{I};
      N = size(SegP,1);
      if SP >= N-3 % Use the whole parent
        Segs{SegInd} = [SegP; SegC];
        Segs{I} = cell(0);
        SPar(I,1) = 0;
        UnMod(I) = false;

        ChildSegs = SChi{I};
        if size(ChildSegs,1) < size(ChildSegs,2)
          ChildSegs = ChildSegs';
        end
        c = set_difference(SChi{SegInd},I,Fal);
        SChi{SegInd} = [c; ChildSegs];
        SPar(ChildSegs,1) = SegInd;
        SPar(ChildSegs,2) = N+SPar(ChildSegs,2);

      else % divide the parent segment into two parts
        ns = ns+1;
        Segs{ns} = SegP(SP+1:end);
        SPar(ns,:) = [SegInd SP];
        Segs{SegInd} = [SegP(1:SP); SegC];
        Segs{I} = cell(0);
        SPar(I,1) = 0;
        UnMod(ns) = true;
        UnMod(I) = false;

        ChildSegs = SChi{SegInd};
        K = SPar(ChildSegs,2) > SP;
        SChi{SegInd} = [ChildSegs(~K); ns];
        ChildSegs = ChildSegs(K);
        SChi{ns} = ChildSegs;
        SPar(ChildSegs,1) = ns;
        SPar(ChildSegs,2) = SPar(ChildSegs,2)-SP;

        ChildSegs = SChi{I};
        c = set_difference(SChi{SegInd},I,Fal);
        SChi{SegInd} = [c; ChildSegs];
        SPar(ChildSegs,1) = SegInd;
        SPar(ChildSegs,2) = SP+SPar(ChildSegs,2);

      end

    end
    UnMod(SegInd) = false;
  else
    UnMod(SegInd) = false;
  end

  % Select the next branch, use increasing branching order
  if BranchOrder > 0 && any(UnMod(SegChildren))
    ChildSegInd = ChildSegInd+1;
    SegInd = SegChildren(ChildSegInd);
  elseif BranchOrder == 0
    BranchOrder = BranchOrder+1;
    SegChildren = SChi{1};
    if ~isempty(SegChildren)
      SegInd = SegChildren(1);
    else
      UnMod = false;
    end
  else
    BranchOrder = BranchOrder+1;
    i = 1;
    SegChildren = SChi{1};
    while i < BranchOrder && ~isempty(SegChildren)
      i = i+1;
      L = cellfun('length',SChi(SegChildren));
      Keep = L > 0;
      SegChildren = SegChildren(Keep);
      SegChildren = vertcat(SChi{SegChildren});
    end
    I = UnMod(SegChildren);
    if any(I)
      SegChildren = SegChildren(I);
      SegInd = SegChildren(1);
      ChildSegInd = 1;
    end
  end
end

% Modify indexes by removing empty segments
Empty = true(ns,1);
for i = 1:ns
  if isempty(Segs{i})
    Empty(i) = false;
  end
end
Segs = Segs(Empty);
Ind = (1:1:ns)';
n = nnz(Empty);
I = (1:1:n)';
Ind(Empty) = I;
SPar = SPar(Empty,:);
J = SPar(:,1) > 0;
SPar(J,1) = Ind(SPar(J,1));
for i = 1:ns
  if Empty(i)
    ChildSegs = SChi{i};
    if ~isempty(ChildSegs)
      ChildSegs = Ind(ChildSegs);
      SChi{i} = ChildSegs;
    end
  end
end
SChi = SChi(Empty);
ns = n;

% Modify SChi
for i = 1:ns
  ChildSegs = SChi{i};
  if size(ChildSegs,2) > size(ChildSegs,1) && size(ChildSegs,1) > 0
    SChi{i} = ChildSegs';
  elseif size(ChildSegs,1) == 0 || size(ChildSegs,2) == 0
    SChi{i} = zeros(0,1);
  end
  Seg = Segs{i};
  n = max(size(Seg));
  for j = 1:n
    ChildSegs = Seg{j};
    if size(ChildSegs,2) > size(ChildSegs,1) && size(ChildSegs,1) > 0
      Seg{j} = ChildSegs';
    elseif size(ChildSegs,1) == 0 || size(ChildSegs,2) == 0
      Seg{j} = zeros(0,1);
    end
  end
  Segs{i} = Seg;
end
end % End of function


function [Segs,SPar,SChi] = remove_small(Ce,Segs,SPar,SChi)

% Removes small child segments

% computes and estimate for stem radius at the base
Segment = Segs{1};  % current or parent segment
ns = size(Segment,1);  % layers in the parent
if ns > 10
  EndL = 10;  % ending layer index in parent
else
  EndL = ns;
end
End = average(Ce(Segment{EndL},:)); % Center of end layer
Start = average(Ce(Segment{1},:));  % Center of starting layer
V = End-Start;  % Vector between starting and ending centers
V = V/norm(V);  % normalize
Sets = vertcat(Segment{1:EndL});
MaxRad = max(distances_to_line(Ce(Sets,:),V,Start));

Nseg = size(Segs,1);
Fal = false(Nseg,1);
Keep = true(Nseg,1);
Sets = zeros(2000,1);
for i = 1:Nseg
  if Keep(i)
    ChildSegs = SChi{i};  % child segments
    if ~isempty(ChildSegs) % child segments exists
      n = length(ChildSegs); % number of children
      Segment = Segs{i};  % current or parent segment
      ns = size(Segment,1);  % layers in the parent
      for j = 1:n % check each child separately
        nl = SPar(ChildSegs(j),2);  % the index of the layer in the parent the child begins
        if nl > 10
          StartL = nl-10; % starting layer index in parent
        else
          StartL = 1;
        end
        if ns-nl > 10
          EndL = nl+10;  % end layer index in parent
        else
          EndL = ns;
        end
        End = average(Ce(Segment{EndL},:));
        Start = average(Ce(Segment{StartL},:));
        V = End-Start;  % Vector between starting and ending centers
        V = V/norm(V);  % normalize

        % cover sets of the child
        ChildSets = Segs{ChildSegs(j)};
        NL = size(ChildSets,1);
        a = 1;
        for k = 1:NL
          S = ChildSets{k};
          Sets(a:a+length(S)-1) = S;
          a = a+length(S);
        end
        ChildSets = Sets(1:a-1);

        % maximum distance in child
        distChild = max(distances_to_line(Ce(ChildSets,:),V,Start));

        if distChild < MaxRad+0.06

          % Select the cover sets of the parent between centers
          NL = EndL-StartL+1;
          a = 1;
          for k = 1:NL
            S = Segment{StartL+(k-1)};
            Sets(a:a+length(S)-1) = S;
            a = a+length(S);
          end
          ParentSets = Sets(1:a-1);

          % maximum distance in parent
          distPar = max(distances_to_line(Ce(ParentSets,:),V,Start));
          if (distChild-distPar < 0.02) || (distChild/distPar < 1.2 && distChild-distPar < 0.06)
            ChildChildSegs = SChi{ChildSegs(j)};
            nc = length(ChildChildSegs);
            if nc == 0
              % Remove, no child segments
              Keep(ChildSegs(j)) = false;
              Segs{ChildSegs(j)} = zeros(0,1);
              SPar(ChildSegs(j),:) = zeros(1,2);
              SChi{i} = set_difference(ChildSegs,ChildSegs(j),Fal);
            else
              L = SChi(ChildChildSegs);
              L = vertcat(L{:}); % child child segments
              if isempty(L)
                J = false(nc,1);
                for k = 1:nc
                  segment = Segs{ChildChildSegs(k)};
                  if isempty(segment)
                    J(k) = true;
                  else
                    segment1 = [vertcat(segment{:}); ParentSets];
                    distSeg = max(distances_to_line(Ce(segment1,:),V,Start));
                    if (distSeg-distPar < 0.02) || (distSeg/distPar < 1.2 && distSeg-distPar < 0.06)
                      J(k) = true;
                    end
                  end
                end
                if all(J)
                  % Remove
                  ChildChildSegs1 = [ChildChildSegs; ChildSegs(j)];
                  nc = length(ChildChildSegs1);
                  Segs(ChildChildSegs1) = cell(nc,1);
                  Keep(ChildChildSegs1) = false;
                  SPar(ChildChildSegs1,:) = zeros(nc,2);
                  d = set_difference(ChildSegs,ChildSegs(j),Fal);
                  SChi{i} = d;
                  SChi(ChildChildSegs1) = cell(nc,1);
                end
              end
            end
          end
        end
      end
    end
    if i == 1
      MaxRad = MaxRad/2;
    end
  end
end
% Modify segments and their indexing
Segs = Segs(Keep);
n = nnz(Keep);
Ind = (1:1:Nseg)';
J = (1:1:n)';
Ind(Keep) = J;
Ind(~Keep) = 0;
SPar = SPar(Keep,:);
J = SPar(:,1) > 0;
SPar(J,1) = Ind(SPar(J,1));
% Modify SChi
for i = 1:Nseg
  if Keep(i)
    ChildSegs = SChi{i};
    if ~isempty(ChildSegs)
      ChildSegs = nonzeros(Ind(ChildSegs));
      if size(ChildSegs,1) < size(ChildSegs,2)
        SChi{i} = ChildSegs';
      else
        SChi{i} = ChildSegs;
      end
    else
      SChi{i} = zeros(0,1);
    end
  end
end
SChi = SChi(Keep);
end % End of function


function [SegP,Base] = modify_parent(P,Bal,Ce,SegP,SegC,nl,PatchDiam,base)

% Expands the base of the branch backwards into its parent segment and
% then removes the expansion from the parent segment.

Base = SegC{1};
if ~isempty(Base)

  % Define the directions of the segments
  DirChi = segment_direction(Ce,SegC,1);
  DirPar = segment_direction(Ce,SegP,nl);

  if length(Base) > 1
    BaseCent = average(Ce(Base,:));
    db = distances_to_line(Ce(Base,:), DirChi', BaseCent); % distances of the sets in the base to the axis of the branch
    DiamBase = 2*max(db);  % diameter of the base
  elseif length(Bal{Base}) > 1
    BaseCent = average(P(Bal{Base},:));
    db = distances_to_line(P(Bal{Base},:), DirChi', BaseCent);
    DiamBase = 2*max(db);
  else
    BaseCent = Ce(Base,:);
    DiamBase = 0;
  end

  % Determine the number of cover set layers "n" to be checked
  Angle = abs(DirChi'*DirPar);  % abs of cosine of the angle between component and segment directions
  Nlayer = max([3,ceil(Angle*2*DiamBase/PatchDiam)]);
  if Nlayer > nl  % can go only to the bottom of the segment
    Nlayer = nl;
  end

  % Check the layers
  layer = 0;
  base{1} = Base;
  while layer < Nlayer
    Sets = SegP{nl-layer};
    Seg = average(Ce(Sets,:)); % mean of the cover sets' centers

    VBase = mat_vec_subtraction(Ce(Sets,:),BaseCent);  % vectors from base's center to sets in the segment
    h = VBase*DirChi;
    B = repmat(DirChi',length(Sets),1);
    B = [h.*B(:,1) h.*B(:,2) h.*B(:,3)];
    V = VBase-B;
    distSets = sqrt(sum(V.*V,2)); % distances of the sets in the segment to the axis of the branch

    VSeg = mat_vec_subtraction(Ce(Sets,:),Seg);  % vectors from segments's center to sets in the segment
    lenBase = sqrt(sum(VBase.*VBase,2)); % lengths of VBase
    lenSeg = sqrt(sum(VSeg.*VSeg,2)); % lengths of VSeg
    if Angle < 0.9
      K = lenBase < 1.1/(1-0.5*Angle^2)*lenSeg;     % sets closer to the base's center than segment's center
      J = distSets < 1.25*DiamBase;   % sets close enough to the axis of the branch
      I = K&J;
    else % branch almost parallel to parent
      I = distSets < 1.25*DiamBase; % only the distance to the branch axis counts
    end

    if all(I) || ~any(I) % stop the process if all the segment's or no segment's sets
      layer = Nlayer;
    else
      SegP{nl-layer} = Sets(not(I));
      base{layer+2} = Sets(I);
      layer = layer+1;
    end
  end
  Base = vertcat(base{1:Nlayer+1});
end

end % End of function


function D = segment_direction(Ce,Seg,nl)

% Defines the direction of the segment

% Define bottom and top layers
if nl-3 > 0
  bot = nl-3;
else
  bot = 1;
end
j = 1;
while j < 3 && isempty(Seg{bot})
  bot = bot+1;
  j = j+1;
end
if nl+2 <= size(Seg,1)
  top = nl+2;
else
  top = size(Seg,1);
end
j = 1;
while j < 3 && isempty(Seg{top})
  top = top-1;
  j = j+1;
end

% Direction
if top > bot
  Bot = average(Ce(Seg{bot},:));
  Top = average(Ce(Seg{top},:));
  V = Top-Bot;
  D = V'/norm(V);
else
  D = zeros(3,1);
end


end % End of function
