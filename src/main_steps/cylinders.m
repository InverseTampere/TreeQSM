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
% Version 2.00
% Latest update     16 Aug 2017
%
% Copyright (C) 2013-2017 Pasi Raumonen
% ---------------------------------------------------------------------

% Reconstructs the surface and volume of branches of input tree with 
% cylinders. Subdivides each segment to smaller regions to which cylinders 
% are fitted in least squares sense. Returns the cylinder information and 
% in addition the child-relation of the cylinders plus the cylinders in 
% each segment.
%
% Inputs:
% P         Point cloud, matrix
% cover     Cover sets
% segment   Segments
% input     Input parameters of the reconstruction, particularly the
%               following
%   lcyl      Cylinder length/radius ratio
%   FilRad    Filtering radius relative to the estimated radius of each
%               subregion. Points further than this are removed from the
%               region and thus from the cylinder fitting.
%
% Outputs:
% cylinder  Structure array containing the following cylinder info: 
%   radius (Rad)        Radii of the cylinders, vector
%   length (Len)        Lengths of the cylinders, vector
%   axis (Axe)          Axes of the cylinders, matrix
%   start (Sta)         Starting points of the cylinders, matrix
%   parent (CPar)       Parents of the cylinders, vector
%   extension (CExt)    Extensions of the cylinders, vector
%   added (Added)       Added cylinders, logical vector
%   UnModRadius (Rad0)  Unmodified radii
%   CylsInSegment       Cylinders (indexes) in each branch, cell-array
%   ChildCyls           Child cylinders of each cylinder, cell-array

Segs = segment.segments;
SPar = segment.ParentSegment;
SChi = segment.ChildSegment;

%% Initialization of variables
NumOfSeg = max(size(Segs));   % number of segments
n = min(40*NumOfSeg,2e5);
Rad = zeros(n,1); % Radii of the cylinders
Len = zeros(n,1); % Lengths of the cylinders
Axe = zeros(n,3); % Axes of the cylinders
Sta = zeros(n,3); % Staring points of the cylinders
CPar = zeros(n,1); % Parents of the cylinders
CExt = zeros(n,1); % Extensions of the cylinders
CChi = cell(n,1); % Children of the cylinders
CiS = cell(NumOfSeg,1); % Cylinders in the segment
SoC = zeros(n,1); % Segment of the cylinder
Added = false(n,1); % is the cylinder added to fill a gap between the branch and parent
Rad0 = zeros(n,1); % Unmodified radii of the cylinders
c = 1;  % number of cylinders determined
Fal = false(size(P,1),1); % auxiliary variable

%% Determine suitable order of segments (from trunk to the "youngest" child)
S = 1;
nc = 1;
SegmentIndex = zeros(NumOfSeg,1);
SegmentIndex(1) = 1;
while ~isempty(S)
    S = vertcat(SChi{S});
    n = length(S);
    SegmentIndex(nc+1:nc+n) = S;
    nc = nc+n;
end

%% Fit cylinders individually for each segment
for k = 1:NumOfSeg
    si = SegmentIndex(k);
    Pass = true;
    if k > 1
        if isempty(CiS{SPar(si)}) || isempty(Segs{si}) || si == 0
            Pass = false;
        end
    end
    if Pass
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
            %% Determine the regions for cylinder fitting
            [Regs,cyl] = regions(P,Points,IndPoints,Seg,inputs,Fal);
            nc = numel(cyl.rad);
            
            %% Fit cylinders to the regions
            if nc > 0
                cyl = cylinder_fitting(P,Regs,cyl);
                nc = numel(cyl.rad);
            end      
            
            %% Search possible parent cylinder
            if nc > 0 && si > 1
                [PC,cyl,added] = parent_cylinder(SPar,SChi,CiS,Rad,Len,Sta,Axe,cyl,si);
                nc = numel(cyl.rad);
            elseif si == 1
                PC = zeros(0,1);
                added = false;
            else
                added = false;
            end
            
            %% Adjust cylinders
            if nc > 0
                cyl = adjustments(Rad,Len,Sta,Axe,cyl,PC,si,inputs);
            end
            
            %% Save the cylinders
            % if at least one acceptable cylinder, then save them
            I = sum(cyl.axe.*cyl.axe,2);
            J = sum(cyl.sta.*cyl.sta,2);
            Accept = nc > 0 && min(cyl.rad(1:nc)) > 0 && ~any(I == 0) && ~any(J == 0);
            if Accept
                % If the parent cylinder exists, set the parent-child relations
                if ~isempty(PC)
                    CPar(c) = PC;
                    if CExt(PC) == c
                        I = SoC(PC);
                        SoC(c:c+nc-1) = I;
                        CiS{I} = [CiS{I}; linspace(c,c+nc-1,nc)'];
                    else
                        CChi{PC} = [CChi{PC}; c];
                        SoC(c:c+nc-1) = si;
                        CiS{si} = linspace(c,c+nc-1,nc)';
                    end
                else
                    SoC(c:c+nc-1) = si;
                    CiS{si} = linspace(c,c+nc-1,nc)';
                end
                
                Rad(c:c+nc-1) = cyl.rad(1:nc);
                Len(c:c+nc-1) = cyl.len(1:nc);
                Axe(c:c+nc-1,:) = cyl.axe(1:nc,:);
                Sta(c:c+nc-1,:) = cyl.sta(1:nc,:);
                CPar(c+1:c+nc-1) = linspace(c,c+nc-2,nc-1);
                CExt(c:c+nc-2) = linspace(c+1,c+nc-1,nc-1);
                Rad0(c:c+nc-1) = cyl.rad0(1:nc);
                if added
                    Added(c) = true;
                end          
                c = c+nc; % number of cylinders so far (plus one)
                
            end
        end
    end
end
c = c-1; % number of cylinders 

%% Define outputs
clear cylinder
cylinder.radius = single(Rad(1:c));     
cylinder.length = single(Len(1:c)); 
cylinder.start = single(Sta(1:c,:));  
cylinder.axis = single(Axe(1:c,:));
if c <= 2^16
    cylinder.parent = uint16(CPar(1:c));
    cylinder.extension = uint16(CExt(1:c));
else
    cylinder.parent = uint32(CPar(1:c));
    cylinder.extension = uint32(CExt(1:c));
end
cylinder.added = logical(Added(1:c));
cylinder.UnmodRadius = single(Rad0(1:c));
for si = 1:NumOfSeg
    if size(CiS{si},2) > 1
        CiS{si} = CiS{si}';
    end
end
cylinder.CylsInSegment = CiS;
CChi = CChi(1:c,:);
cylinder.ChildCyls = CChi;

% Growth volume correction
if inputs.GrowthVolCor && c > 0
    cylinder = growth_volume_correction(cylinder,inputs);
end

end % End of main function


function [Regs,cyl] = regions(P,Points,Ind,Seg,inputs,Fal)

% Define the subregions of the segment for cylinder fitting

nl = max(size(Seg)); % number of cover set layers in the segment
if nl > 3
    %% Define each region with approximate relative length of lcyl

    % Define first region
    Fal(Points) = true;
    Test = Points(Ind(1,1):Ind(4,2));
    Bot = Points(Ind(1,1):Ind(2,2));
    Bot = average(P(Bot,:));
    Top = Points(Ind(3,1):Ind(4,2));
    Top = average(P(Top,:));
    V = Top-Bot;
    [d,~,h] = distances_to_line(P(Test,:),V,Bot);
    [~,~,ht] = distances_to_line(Top,V,Bot);
    I = h <= ht;
    R = median2(d(I));
    J = d < inputs.FilRad*R;
    I = I&J;
    R = median2(d(I));
    L = max(h(I))-min(h);
    i = 4;   i0 = 1;
    while (i <= nl-1 && L < inputs.lcyl*R) || (i <= nl-1 && nnz(I) < 30)
        top = Points(Ind(i,1):Ind(i,2));
        Top = average(P(top,:));
        V = Top-Bot;
        i = i+1;
        Test = Points(Ind(i0,1):Ind(i,2));
        [d,~,h] = distances_to_line(P(Test,:),V,Bot);
        [~,~,ht] = distances_to_line(Top,V,Bot);
        I = h <= ht;
        R = median2(d(I));
        J = d < inputs.FilRad*R;
        I = I&J;
        R = median2(d(I));
        L = max(h(I))-min(h);
    end
    Test = Test(I);
    if i == 4 && L/R > inputs.lcyl
        NL = 3;
    else
        NL = i;
    end
    % Initialize regions and cylinders
    n = ceil(3*nl/NL);
    Regs = cell(n,1);       
    Regs{1} = Test;
    Fal(Test) = false;
    Axes = zeros(n,3);      
    Axes(1,:) = V'/norm(V); 
    Starts = zeros(n,3);    
    Starts(1,:) = Bot+min(h)*Axes(1,:);
    Rads = zeros(n,1);      
    Rads(1) = R;
    Lengs = zeros(n,1);     
    Lengs(1) = L;
    
    % Define the other regions
    t = 1;
    i0 = NL-1;
    i = NL+ceil(NL/3);
    while i <= nl-1
        k = ceil(NL/3);
        Bot = Top;
        top = Points(Ind(i,1):Ind(i,2));
        Top = average(P(top,:));
        V = Top-Bot;
        i = i+1;
        Test = Points(Ind(i0,1):Ind(i,2));
        [d,~,h] = distances_to_line(P(Test,:),V,Bot);
        [~,~,ht] = distances_to_line(Top,V,Bot);
        I = h <= ht;
        J = h >= 0;
        I = I&J;
        if nnz(I) < 3
            I = h >= 0;
        end
        R = median2(d(I));
        if R == 0
            R = average(d(I));
            if R == 0
                R = max(d(I));
                if R == 0
                    R = mad(d);
                end
            end
        end
        J = d < inputs.FilRad*R;
        I = I&J;
        R = median2(d(I));
        L = norm(V);
        k = k+1;
        while (i <= nl-1 && L < inputs.lcyl*R && k <= NL) || (i <= nl-1 && nnz(I) < 20)
            top = Points(Ind(i,1):Ind(i,2));
            Top = average(P(top,:));
            V = Top-Bot;
            i = i+1;
            Test = Points(Ind(i0,1):Ind(i,2));
            [d,~,h] = distances_to_line(P(Test,:),V,Bot);
            [~,~,ht] = distances_to_line(Top,V,Bot);
            I = h <= ht;
            J = h >= 0;
            I = I&J;
            if nnz(I) < 3
                I = h >= 0;
            end
            R = median2(d(I));
            if R == 0
                R = mean(d(I));
                if R == 0
                    R = max(d(I));
                    if R == 0
                        R = mad(d);
                    end
                end
            end
            J = d < inputs.FilRad*R;
            I = I&J;
            R = median2(d(I));
            L = norm(V);
            k = k+1;
        end
        if i >= nl-1
            Test = Points(Ind(i0,1):Ind(nl,2));
            [d,~,h] = distances_to_line(P(Test,:),V,Bot);
            I = h >= 0;
            R = median2(d(I));
            if R == 0
                R = mean(d(I));
                if R == 0
                    R = max(d(I));
                    if R == 0
                        R = mad(d);
                    end
                end
            end
            J = d < inputs.FilRad*R;
            I = I&J;
            R = median2(d(I));
            L = max(h(I));
            Test = Test(I);
            I = Fal(Test);
            Test = Test(I);
            if length(Test) >= 20
                t = t+1;
                Regs{t} = Test;
                Axes(t,:) = V'/norm(V);     
                Starts(t,:) = Bot;
                Rads(t) = R;                
                Lengs(t) = L;
            else
                Regs{t} = [Regs{t}; Test];
            end
        else
            Test = Test(I);
            I = Fal(Test);
            Test = Test(I);
            if length(Test) >= 20
                t = t+1;
                Regs{t} = Test;
                Axes(t,:) = V'/norm(V);     
                Starts(t,:) = Bot;
                Rads(t) = R;                
                Lengs(t) = L;
            else
                Regs{t} = [Regs{t}; Test];
            end
        end
        i0 = i-1;
        i = i0+ceil(NL/3);
    end
    Axes = Axes(1:t,:);     
    V = Starts(2:t,:)-Starts(1:t-1,:);      
    L = sqrt(sum(V.*V,2));
    Axes(1:t-1,:) = [V(:,1)./L V(:,2)./L V(:,3)./L];
    Starts = Starts(1:t,:);         
    Rads = Rads(1:t);
    Lengs = Lengs(1:t);     
    Lengs(1:t-1) = L;
    Regs = Regs(1:t);
    
else
    %% Define a region for small segments
    % Define the direction
    Bot = Points(Ind(1,1):Ind(1,2));
    Bot = average(P(Bot,:));
    Top = Points(Ind(nl,1):Ind(nl,2));
    Top = average(P(Top,:));
    Axes = Top-Bot;
    Axes = Axes/norm(Axes);
        
    % Define other outputs
    Regs = cell(1,1);
    Regs{1} = Points;
    Starts = average(P(Points,:));
    if max(size(Starts)) == 3
        [d,~,h] = distances_to_line(P(Points,:),Axes,Starts);
        Lengs = max(h)-min(h);
        R = median2(d);
        I = d < inputs.FilRad*R;
        Rads = median2(d(I));
        Height = P(Points,:)*Axes';
        hpoint = Starts*Axes';
        Starts = Starts-(hpoint-min(Height))*Axes;
        t = 1;
    else
        t = 0;
        Axes = 0;
        Rads = 0;
        Lengs = 0;
    end
    
end

if (t > 1) && (length(Regs{t}) < 11)
    Regs{t-1} = [Regs{t-1}; Regs{t}];
    t = t-1;
    Regs = Regs(1:t);
    Rads = Rads(1:t);
    Lengs = Lengs(1:t);
    Axes = Axes(1:t,:);
    Starts = Starts(1:t,:);
end

clear cyl
cyl.rad = Rads;
cyl.len = Lengs;
cyl.sta = Starts;
cyl.axe = Axes;
cyl.rad0 = Rads;
cyl.len0 = Lengs;
cyl.sta0 = Starts;
cyl.axe0 = Axes;

end % End of function


function cyl = cylinder_fitting(P,Regs,cyl)

% Fit cylinders to the regions
warning off
nr = size(Regs,1); % number of regions
ci = 0; % cylinder index
for j = 1:nr
    if (length(Regs{j}) > 10) && (norm(cyl.axe0(j,:)) > 0) % fit cylinders to large enough subsegs
        
        % Initial estimates
        Region = P(Regs{j},:);  % the coordinate points used for fitting
        Axis0 = cyl.axe0(j,:);     % cylinder axis
        Point0 = cyl.sta0(j,:);  % point in the cylinder axis
        R0 = cyl.rad0(j,1);
        
        %% First fitting
        [R,L,Point,Axis,d,conv,rel] = least_squares_cylinder(Region,Point0,Axis0,R0);
        
        % Conditions for second fitting and accepting the results
        I1 = conv & rel; % fitting converged and is reliable
        I2 = ~(isnan(R)|any(isnan(Point))|any(isnan(Axis))); % results are numbers
        mad = average(abs(d));  % mean distance to the fitted cylinder
        md = max(d);  % maximum distance
        I3 = mad < R0 & abs(Axis0*Axis') > 0.8; % distances and the angle acceptable
        I4 = R < 3*R0 & R > 0; % radius is acceptable
        % second fitting if large enough absolute and relative "errors"
        SecondFitting = mad > 0.005 & mad/R > 0.05 & md/R > 0.2;
        AcceptFitting = I1&I2&I3&I4; % accept the first fitting
        SecondFitting = AcceptFitting & SecondFitting; % second cylinder fitting
        
        % Possible second fitting
        if SecondFitting
            
            if mad > 0.015 && R > 0.03
                %% Try two shorter cylinders
                % If large cylinder with large average error/distance, try
                % replacing the cylinder with two shorter ones. Try five
                % different combinations: 3-7, 4-6, 5-5, 6-4, 7-3 portions.
                h = Region*Axis';
                hmin = min(h);
                
                % Best results
                madb = mad;
                Rb = R;  Lb = L;  Pointb = Point;  Axisb = Axis;
                
                % Try all the 5 combinations, save always the best one
                for i = 1:5
                    % Define the two subregions
                    I = h <= hmin+(0.2+0.1*i)*L;
                    Region1 = Region(I,:);
                    Region2 = Region(~I,:);
                    if nnz(I) > 10 && nnz(~I) > 10
                        
                        % Fit cylinders
                        [R1,L1,Point1,Axis1,d1,conv1,rel1] = least_squares_cylinder(Region1,Point,Axis,R);
                        
                        [R2,L2,Point2,Axis2,d2,conv2,rel2] = least_squares_cylinder(Region2,Point,Axis,R);
                        
                        % Check if acceptable and also the best fits so far
                        if conv1 && rel1 && conv2 && rel2
                            % Fit is the best one if both mad1 and mad2 are
                            % smaller than mad (from the one fitted
                            % cylinder) and if their average is smaller
                            % than the best madb
                            mad1 = average(abs(d1));
                            mad2 = average(abs(d2));
                            AcceptFits = mad1 < mad && mad2 < mad && (mad1+mad2)/2 < madb;
                            AcceptFits = AcceptFits && abs(Axis*Axis1') > 0.8 && abs(Axis*Axis2') > 0.7;
                            AcceptFits = AcceptFits && R1 < 1.33*R && R2 < 1.33*R;
                        else
                            AcceptFits = false;
                        end
                        
                        if AcceptFits
                            % Update the best results
                            madb = (mad1+mad2)/2;
                            Rb = [R1; R2];
                            Lb = [L1; L2];
                            Pointb = [Point1; Point2];
                            Axisb = [Axis1; Axis2];
                             
                            % Try second fit with outliers removed
                            % Remove the outliers
                            I = d1 < 0;
                            d1(I) = 0.5*abs(d1(I));
                            [~,J] = sort(d1);
                            I = J(1:ceil(0.7*length(J)));
                            k = 0;
                            while length(I) < 10
                                k = k+1;
                                I = J(1:ceil((0.7+0.1*k)*length(J)));
                            end
                            Region1 = Region1(I,:);
                            
                            I = d2 < 0;
                            d2(I) = 0.5*abs(d2(I));
                            [~,J] = sort(d2);
                            I = J(1:ceil(0.7*length(J)));
                            k = 0;
                            while length(I) < 10
                                k = k+1;
                                I = J(1:ceil((0.7+0.1*k)*length(J)));
                            end
                            Region2 = Region2(I,:);
                            
                            % Fit again
                            [R12,L1,Point1,Axis1,d1,conv1,rel1] = least_squares_cylinder(Region1,Point1,Axis1,R1);
                            
                            [R22,L2,Point2,Axis2,d2,conv2,rel2] = least_squares_cylinder(Region2,Point2,Axis2,R2);
                            
                            % Update the best results if ok fits
                            if conv1 && rel1 && conv2 && rel2 && R12 < 1.1*R1 && R22 < 1.1*R2
                                Rb = [R12; R22];
                                Lb = [L1; L2];
                                Pointb = [Point1; Point2];
                                Axisb = [Axis1; Axis2];
                            end
                        end
                    end
                end
                % Use the best results
                R = Rb;  L = Lb;  Point = Pointb;  Axis = Axisb;
            else
                %% Second fit with outliers removed
                % Save the first fit results
                R1 = R;  L1 = L;  Axis1 = Axis;  Point1 = Point;
                
                % Remove the outliers
                I = d < 0;
                d(I) = 0.5*abs(d(I));
                [~,J] = sort(d);
                I = J(1:ceil(0.7*length(J)));
                k = 0;
                while length(I) < 10
                    k = k+1;
                    I = J(1:ceil((0.7+0.1*k)*length(J)));
                end
                Region = Region(I,:);
                
                % Second fitting
                [R,L,Point,Axis,d,conv,rel] = least_squares_cylinder(Region,Point0,Axis0,R0);
                
                % Conditions for accepting the results
                I1 = conv & rel; % fitting converged and is reliable
                I2 = ~(isnan(R)|any(isnan(Point))|any(isnan(Axis))); % results are numbers
                mad = average(abs(d));  % mean distance to the fitted cylinder
                I3 = mad < R0 & abs(Axis0*Axis') > 0.8; % distances and the angle acceptable
                I4 = R < 3*R0 & R > 0; % radius is acceptable
                AcceptFitting = I1&I2&I3&I4; % accept the first fitting
            end
            
            if ~AcceptFitting
                % if the second fit was bad, use the results from the first fit
                R = R1;  L = L1;  Axis = Axis1;  Point = Point1;
            end
            
        end
        
        if AcceptFitting
            % Save the new fitted values
            if length(R) == 1
                ci = ci+1;
                cyl.rad(ci,1) = R;
                cyl.len(ci,1) = L;
                cyl.sta(ci,:) = Point;
                cyl.axe(ci,:) = Axis;
            else
                ci = ci+1;
                cyl.rad(ci:ci+1,1) = R;
                cyl.len(ci:ci+1,1) = L;
                cyl.sta(ci:ci+1,:) = Point;
                cyl.axe(ci:ci+1,:) = Axis;
                cyl.sta0 = [cyl.sta0(1:ci,:); cyl.sta0(ci:end,:)];
                cyl.sta0(ci+1,:) = cyl.sta0(ci,:)+cyl.len0(ci)/2*cyl.axe0(ci,:);
                cyl.axe0 = [cyl.axe0(1:ci,:); cyl.axe0(ci:end,:)];
                cyl.rad0 = [cyl.rad0(1:ci,:); cyl.rad0(ci:end,:)];
                cyl.len0 = [cyl.len0(1:ci,:)/2; cyl.len0(ci:end,:)/2];
                ci = ci+1;
            end
        else
            % do not accept least square fittings, use initial estimates
            ci = ci+1;
            cyl.rad(ci,1) = cyl.rad0(ci,1);
            cyl.len(ci,1) = cyl.len0(ci,1);
            cyl.sta(ci,:) = cyl.sta0(ci,:);
            cyl.axe(ci,:) = cyl.axe0(ci,:);
        end
    end
end
warning on

% if ci > 0 %nr-3
%     figure(6)
%     subplot(1,2,1)
%     plot_segs(P,Regs,6,5)
%     hold on
%     for i = 1:ci
%       draw(cyl.rad(i),cyl.len(i),cyl.axe(i,:),cyl.sta(i,:),1,20)
%     end
%     hold off
%     subplot(1,2,2)
%     plot_segs(P,Regs,6,5)
%     hold on
%     for i = 1:ci
%       draw(cyl.rad0(i),cyl.len0(i),cyl.axe0(i,:),cyl.sta0(i,:),1,20)
%     end
%     hold off
%     pause
% end

end % End of function


function [PC,cyl,added] = parent_cylinder(SPar,SChi,CiS,Rad,Len,Sta,Axe,cyl,si)

% Finds the parent cylinder from the possible parent segment.
% Does this by checking if the axis of the cylinder, if continued, will
% cross the nearby cylinders in the parent segment.
% Adjust the cylinder so that it starts from the surface of its parent.

Rads = cyl.rad;
Lengs = cyl.len;
Starts = cyl.sta;
Axes = cyl.axe;

% PC     Parent cylinder
nc = numel(Rads);
added = false;
if SPar(si) > 0 % parent segment exists, find the parent cylinder
    s = SPar(si);
    PC = CiS{s}; % the cylinders in the parent segment
    % select the closest cylinders for closer examination
    if length(PC) > 1
        D = mat_vec_subtraction(-Sta(PC,:),-Starts(1,:));
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
        for j = 1:n
            % Crossing points solved from a quadratic equation
            A = Axes(1,:)-(Axes(1,:)*Axe(pc(j),:)')*Axe(pc(j),:);
            B = Starts(1,:)-Sta(pc(j),:)-(Starts(1,:)*Axe(pc(j),:)')*Axe(pc(j),:)...
                +(Sta(pc(j),:)*Axe(pc(j),:)')*Axe(pc(j),:);
            e = A*A';
            f = 2*A*B';
            g = B*B'-Rad(pc(j))^2;
            di = sqrt(f^2 - 4*e*g);  % the discriminant
            s1 = (-f + di)/(2*e);       
            s2 = (-f - di)/(2*e); % how much the starting point must be moved to cross
            if isreal(s1) %% cylinders can cross
                % the heights of the crossing points
                x(j,:) = [s1 s2];
                h(j,1) = Starts(1,:)*Axe(pc(j),:)'+x(j,1)*Axes(1,:)*Axe(pc(j),:)'-...
                    Sta(pc(j),:)*Axe(pc(j),:)';
                h(j,2) = Starts(1,:)*Axe(pc(j),:)'+x(j,2)*Axes(1,:)*Axe(pc(j),:)'-...
                    Sta(pc(j),:)*Axe(pc(j),:)';
            end
        end
        
        %% Extend to crossing point in the (extended) parent
        I = x(:,1) ~= 0; % Select only candidates with crossing points
        pc = pc0(I);    x = x(I,:);     h = h(I,:);
        j = 1;      n = nnz(I);         
        X = zeros(n,3); % 
        while j <= n && ~ParentFound
            if x(j,1) > 0 && x(j,2) < 0
                % sp inside the parent and crosses its surface
                if h(j,1) >= 0 && h(j,1) <= Len(pc(j)) && Lengs(1)-x(j,1) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,1)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,1);
                    ParentFound = true;
                elseif Lengs(1)-x(j,1) > 0
                    if h(j,1) < 0
                        X(j,:) = [x(j,1) abs(h(j,1)) 0];
                    else
                        X(j,:) = [x(j,1) h(j,1)-Len(pc(j)) 0];
                    end
                else
                    X(j,:) = [x(j,1) h(j,1) 1];
                end
            elseif x(j,1) < 0 && x(j,2) > 0 && Lengs(1)-x(j,2) > 0
                % sp inside the parent and crosses its surface
                if h(j,2) >= 0 && h(j,2) <= Len(pc(j)) && Lengs(1)-x(j,2) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,2)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,2);
                    ParentFound = true;
                elseif Lengs(1)-x(j,2) > 0
                    if h(j,2) < 0
                        X(j,:) = [x(j,2) abs(h(j,2)) 0];
                    else
                        X(j,:) = [x(j,2) h(j,2)-Len(pc(j)) 0];
                    end
                else
                    X(j,:) = [x(j,2) h(j,2) 1];
                end
            elseif x(j,1) < 0 && x(j,2) < 0 && x(j,2) < x(j,1) && Lengs(1)-x(j,1) > 0
                % sp outside the parent and crosses its surface when extended
                % backwards
                if h(j,1) >= 0 && h(j,1) <= Len(pc(j)) && Lengs(1)-x(j,1) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,1)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,1);
                    ParentFound = true;
                elseif Lengs(1)-x(j,1) > 0
                    if h(j,1) < 0
                        X(j,:) = [x(j,1) abs(h(j,1)) 0];
                    else
                        X(j,:) = [x(j,1) h(j,1)-Len(pc(j)) 0];
                    end
                else
                    X(j,:) = [x(j,1) h(j,1) 1];
                end
            elseif x(j,1) < 0 && x(j,2) < 0 && x(j,2) > x(j,1) && Lengs(1)-x(j,2) > 0
                % sp outside the parent and crosses its surface when extended
                % backwards
                if h(j,2) >= 0 && h(j,2) <= Len(pc(j)) && Lengs(1)-x(j,2) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,2)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,2);
                    ParentFound = true;
                elseif Lengs(1)-x(j,2) > 0
                    if h(j,2) < 0
                        X(j,:) = [x(j,2) abs(h(j,2)) 0];
                    else
                        X(j,:) = [x(j,2) h(j,2)-Len(pc(j)) 0];
                    end
                else
                    X(j,:) = [x(j,2) h(j,2) 1];
                end
            elseif x(j,1) > 0 && x(j,2) > 0 && x(j,2) < x(j,1) && Lengs(1)-x(j,1) > 0
                % sp outside the parent but crosses its surface when extended
                % forward
                if h(j,1) >= 0 && h(j,1) <= Len(pc(j)) && Lengs(1)-x(j,1) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,1)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,1);
                    ParentFound = true;
                elseif Lengs(1)-x(j,1) > 0
                    if h(j,1) < 0
                        X(j,:) = [x(j,1) abs(h(j,1)) 0];
                    else
                        X(j,:) = [x(j,1) h(j,1)-Len(pc(j)) 0];
                    end
                else
                    X(j,:) = [x(j,1) h(j,1) 1];
                end
            elseif x(j,1) > 0 && x(j,2) > 0 && x(j,2) > x(j,1) && Lengs(1)-x(j,2) > 0
                % sp outside the parent and crosses its surface when extended
                % forward
                if h(j,2) >= 0 && h(j,2) <= Len(pc(j)) && Lengs(1)-x(j,2) > 0
                    PC = pc(j);
                    Starts(1,:) = Starts(1,:)+x(j,2)*Axes(1,:);
                    Lengs(1) = Lengs(1)-x(j,2);
                    ParentFound = true;
                elseif Lengs(1)-x(j,2) > 0
                    if h(j,1) < 0
                        X(j,:) = [x(j,2) abs(h(j,2)) 0];
                    else
                        X(j,:) = [x(j,2) h(j,2)-Len(pc(j)) 0];
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
            if X(3) == 0 && H < 0.1*Len(pc(I))
                PC = pc(I);
                Starts(1,:) = Starts(1,:)+X(1)*Axes(1,:);
                Lengs(1) = Lengs(1)-X(1);
                ParentFound = true;
            else
                PC = pc(I);
                
                if nc > 1 && X(1) <= Rads(1) && abs(X(2)) <= 1.25*Len(PC)
                    % Remove the first cylinder and adjust the second
                    S = Starts(1,:)+X(1)*Axes(1,:);
                    V = Starts(2,:)+Lengs(2)*Axes(2,:)-S;
                    Lengs(2) = norm(V);         Lengs = Lengs(2:nc);
                    Axes(2,:) = V/norm(V);      Axes = Axes(2:nc,:);
                    Starts(2,:) = S;            Starts = Starts(2:nc,:);
                    Rads = Rads(2:nc);
                    nc = nc-1;
                    ParentFound = true;
                elseif nc > 1
                    % Remove the first cylinder
                    Starts = Starts(2:nc,:);     Lengs = Lengs(2:nc);
                    Axes = Axes(2:nc,:);         Rads = Rads(2:nc);
                    nc = nc-1;
                elseif isempty(SChi{si})
                    % Remove the cylinder
                    nc = 0;
                    PC = zeros(0,1);
                    ParentFound = true;
                    Rads = zeros(0,1);
                elseif X(1) <= Rads(1) && abs(X(2)) <= 1.5*Len(PC)
                    % Adjust the cylinder
                    Starts(1,:) = Starts(1,:)+X(1)*Axes(1,:);
                    Lengs(1) = abs(X(1));
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
            
            [Dist,~,DistOnLines] = distances_between_lines(Starts(1,:),Axes(1,:),Sta(pc,:),Axe(pc,:));
            
            I = DistOnLines >= 0;
            J = DistOnLines <= Len(pc);
            I = I&J;
            if ~any(I)
                I = DistOnLines >= -0.2*Len(pc);
                J = DistOnLines <= 1.2*Len(pc);
                I = I&J;
            end
            if any(I)
                pc = pc(I);     Dist = Dist(I);     DistOnLines = DistOnLines(I);
                [~,I] = min(Dist);
                DistOnLines = DistOnLines(I);       PC = pc(I);
                Q = Sta(PC,:)+DistOnLines*Axe(PC,:);
                V = Starts(1,:)-Q;      L = norm(V);        V = V/L;
                a = acos(V*Axe(PC,:)');
                h = sin(a)*L;
                S = Q+Rad(PC)/h*L*V;
                L = (h-Rad(PC))/h*L;
                if L > 0.01 && L/Lengs(1) > 0.2
                    nc = nc+1;
                    Starts = [S; Starts];       Rads = [Rads(1); Rads];
                    Axes = [V; Axes];       Lengs = [L; Lengs];
                    added = true;
                end
            else
                V = -mat_vec_subtraction(Sta(pc,:),Starts(1,:));
                L = sqrt(sum(V.*V,2));
                V = [V(:,1)./L V(:,2)./L V(:,3)./L];
                A = V*Axes(1,:)';
                [A,I] = max(A);
                L = L(I);       PC = pc(I);     V = V(I,:);
                a = acos(V*Axe(PC,:)');
                h = sin(a)*L;
                S = Sta(PC,:)+Rad(PC)/h*L*V;
                L = (h-Rad(PC))/h*L;
                if L > 0.01 && L/Lengs(1) > 0.2
                    nc = nc+1;
                    Starts = [S; Starts];       Rads = [Rads(1); Rads];
                    Axes = [V; Axes];       Lengs = [L; Lengs];
                    added = true;
                end
            end
        end
    end
else
    % no parent segment exists
    PC = zeros(0,1);
    display('No parent segment')
end

Rads = Rads(1:nc);       
Lengs = Lengs(1:nc,:);
Axes = Axes(1:nc,:);     
Starts = Starts(1:nc,:);

cyl.rad0 = Rads;
if nc > 0 && Lengs(1) ~= cyl.len0(1) && ~added
    cyl.sta0(1,:) = Starts(1,:);
elseif added
    cyl.sta0 = [Starts(1,:); cyl.sta0];
    cyl.axe0 = [Axes(1,:); cyl.axe0];
end

cyl.rad = Rads;
cyl.len = Lengs;
cyl.sta = Starts;
cyl.axe = Axes;

end % End of function


function cyl = adjustments(Rad,Len,Sta,Axe,cyl,PC,si,inputs)

Rads = cyl.rad;
Lengs = cyl.len;
Starts = cyl.sta;
Axes = cyl.axe;
Starts0 = cyl.sta0;
Axes0 = cyl.axe0;

nc = size(Rads,1);
Rads0 = Rads;

MinR = inputs.MinCylRad;
%% Determine the maximum radius based on parent branch
if ~isempty(PC)
    MaxR = 0.95*Rad(PC);
    MaxR = max(MaxR,MinR);
elseif si == 1
    % For the trunk use the maximum from the bottom cylinders
    a = min(3,nc);
    MaxR = 1.25*max(Rads(1:a));
else
    MaxR = 0.005;
end

%% Check maximum and minimum radii
I = Rads < MinR;
Rads(I) = MinR;
if inputs.ParentCor
    I = Rads > MaxR;
    Rads(I) = MaxR;
end

%% Use taper correction to modify radius of too small and large cylinders
if inputs.TaperCor
    if max(Rads) < 0.005
        
        %% Adjust radii of thin branches to be linearly decreasing
        if nc > 2
            r = sort(Rads);
            r = r(2:end-1);
            a = 2*mean(r);
            if a > max(r)
                a = min(0.01,max(r));
            end
            b = min(0.5*min(Rads),0.001);
            Rads = linspace(a,b,nc)';
        else
            r = max(Rads);
            if nc == 1
                Rads = r;
            else
                Rads = [r; 0.5*r];
            end
        end
        
    elseif nc > 4
        %% Parabola adjustment of maximum and minimum
        % Define parabola taper shape as maximum radii
        % "a" is the number first radii used to determine base radius
        r0 = MinR;
        l = sum(Lengs(1:nc)); % branch length
        L = zeros(nc,1); % middle points of cylinder as cumulative length from base
        for i = 1:nc
            if i > 1
                L(i) = Lengs(i)/2+sum(Lengs(1:i-1));
            else
                L(i) = Lengs(i)/2;
            end
        end
        a = 1;
        while L(a) < 0.1*L(end)
            a = a+1;
        end
        a = max(a,2);
        r = 1.05*sum(Rads(1:a))/a; % branch base radius
        
        if si > 1
            % Determine data "S" for parabola fitting
            b = round(nc/4);
            if b >= 3
                %  use 3 first 1/4-length sections to define data points as mean
                %  radii of those cylinders
                S = zeros(5,2);
                S(1,2) = r;
                I0 = 1;
                for i = 1:3
                    [~,I] = min(abs(L-i/4*l));
                    S(i+1,1) = L(I);
                    S(i+1,2) = 1.05*mean(Rads(I0:I));
                    I0 = I+1;
                end
            else
                j = 1;
                S = zeros(5,2);
                S(1,2) = r;
                for i = 1:3
                    S(i+1,1) = L(round(j+b/2));
                    S(i+1,2) = 1.05*sum(Rads(j:j+b-1))/b;
                    j = j+b;
                end
            end
            S(5,:) = [l r0];
            
            % Least square fitting of parabola to "S"
            A = [sum(S(:,1).^4) sum(S(:,1).^2); sum(S(:,1).^2) 5];
            y = [sum(S(:,2).*S(:,1).^2); sum(S(:,2))];
            x = A\y;
            R = x(1)*L.^2+x(2); % parabola
            I = Rads > R;
            Rads(I) = R(I);  % change values larger than parabola-values
            Q = 0.75*R;
            I = Q < r0;
            Q(I) = r0;
            I = Rads < Q;
            Rads(I) = Q(I);
        elseif si == 1 && nc > 5
            % Define partially linear maximum taper curve data S with 8
            % sections
            a = 1;
            while L(a) < 0.06*L(end)
                a = a+1;
            end
            a = max(a,2);
            r = 1.2*max(Rads(1:a));
            b = round(nc/8);
            if b >= 3
                %  use 6 first 1/8-length parts to define data points as mean
                %  radii of those cylinders
                S = zeros(7,2);
                S(1,2) = r;
                I0 = 1;
                for i = 1:5
                    [~,I] = min(abs(L-i/8*l));
                    S(i+1,1) = L(I);
                    S(i+1,2) = 1.05*mean(Rads(I0:I));
                    if S(i+1,2) > S(i,2)
                        S(i+1,2) = S(i,2);
                    end
                    I0 = I+1;
                end
            else
                j = 1;
                S = zeros(7,2);
                S(1,2) = r;
                for i = 1:5
                    S(i+1,1) = L(round(j+b/2));
                    S(i+1,2) = 1.05*sum(Rads(j:j+b-1))/b;
                    j = j+b;
                end
            end
            S(7,:) = [l r0];
            % Check the radii against the taper data (minimum allowed is 70% of maximum)
            j = 1;
            for i = 1:nc
                R = S(j,2)+(L(i)-S(j,1))/(S(j+1,1)-S(j,1))*(S(j+1,2)-S(j,2));
                if Rads(i) > R
                    Rads(i) = R;
                elseif Rads(i) < 0.7*R
                    Rads(i) = 0.7*R;
                end
                if j < 6 && L(i) >= S(j+1,1)
                    j = j+1;
                end
            end
        else
            % Define partially linear maximum taper curve data S with 2
            % sections
            r = 1.1*mean(Rads(1));
            S = zeros(3,2);
            S(1,2) = r;
            S(2,:) = [L(3) 1.05*(Rads(2)+Rads(3))/2];
            S(3,:) = [l r0];
            % Check the radii against the taper data (minimum allowed is 70% of maximum)
            j = 1;
            for i = 1:nc
                R = S(j,2)+(L(i)-S(j,1))*(S(j+1,2)-S(j,2))/(S(j+1,1)-S(j,1));
                if Rads(i) > R
                    Rads(i) = R;
                elseif Rads(i) < 0.7*R
                    Rads(i) = 0.7*R;
                end
                if j < 3 && L(i) > S(j+1,1)
                    j = j+1;
                end
            end
        end
    else
        %% Adjust radii of short branches to be linearly decreasing
        if nc > 2
            a = 2*(Rads(1)+Rads(2))/2;
            if a > max(Rads)
                a = max(Rads);
            end
            b = MinR;
            Rads = linspace(a,b,nc)';
        else
            r = max(Rads);
            if nc == 1
                Rads = r;
            else
                Rads = [r; 0.5*r];
            end
        end
        
    end
end

%% Check big adjustments of starting points
% If modification of the radius was large, then adjust the starting point
% to the initial value
for i = 1:nc
    d = abs(Rads(i)-Rads0(i));
    if d > 0.01 || d > 0.5*Rads0(i)
        S0 = Starts0(i,:);
        S = Starts(i,:);
        A0 = Axes0(i,:);
        d = distances_to_line(S,A0,S0);
        if d > 0.01 || d > 0.5*Rads0(i)
            Starts(i,:) = S0;
            Axes(i,:) = A0;
        end
    end
end

%% Continuous branches
% Make cylinders properly "continuous" by moving the starting points
% First check, move the starting point to the plane defined by parent
% cylinder's top
if nc > 1
    for j = 2:nc
        U = Starts(j,:)-Starts(j-1,:)-Lengs(j-1)*Axes(j-1,:);
        if (norm(U) > 0.0001)
            % First define vector V and W which are orthogonal to the
            % cylinder axis N
            N = Axes(j,:)';
            if norm(N) > 0
                [V,W] = orthonormal_vectors(N);
                % Now define the new starting point
                x = [N V W]\U';
                Starts(j,:) = Starts(j,:)-x(1)*N';
                if x(1) > 0
                    Lengs(j) = Lengs(j)+x(1);
                elseif Lengs(j)+x(1) > 0
                    Lengs(j) = Lengs(j)+x(1);
                end
            end
        end
    end
end

%% Connect far away first cylinders to the parent
if si > 1
    [d,V,h,B] = distances_to_line(Starts(1,:),Axe(PC,:),Sta(PC,:));
    d = d-Rad(PC);
    if d > 0.01
        S = Starts(1,:);
        E = S+Lengs(1)*Axes(1,:);
        V = Rad(PC)*V/norm(V);
        if h >= 0 && h <= Len(PC)
            Starts(1,:) = Sta(PC,:)+B+V;
        elseif h < 0
            Starts(1,:) = Sta(PC,:)+V;
        else
            Starts(1,:) = Sta(PC,:)+Len(PC)*Axe(PC,:)+V;
        end
        Axes(1,:) = E-Starts(1,:);
        Lengs(1) = norm(Axes(1,:));
        Axes(1,:) = Axes(1,:)/Lengs(1);
    end
end

cyl.rad = Rads;
cyl.len = Lengs;
cyl.sta = Starts;
cyl.axe = Axes;

end % End of function
