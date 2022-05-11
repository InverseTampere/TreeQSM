function cylinder = growth_volume_correction(cylinder,inputs)

% ---------------------------------------------------------------------
% GROWTH_VOLUME_CORRECTION.M       Use growth volume allometry approach to 
%                                   modify the radius of cylinders.
%
% Version 2.0.0
% Latest update     16 Sep 2021
%
% Copyright (C) 2013-2021 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Use growth volume (= the total volume "supported by the cylinder") 
% allometry approach to modify the radius of too large and too small 
% cylinders. Uses the allometry: 
%
%       Radius = a * GrowthVolume^b + c
%
% If cylinder's radius is over fac-times or under 1/fac-times the radius
% predicted from the growth volume allometry, then correct the radius to 
% match the allometry. However, the radius of the cylinders in the branch
% tips are never incresed, only decreased by the correction. More details 
% can be from Jan Hackenberg's "SimpleTree" papers and documents.
% ---------------------------------------------------------------------
% Inputs:
% cylinder    Structure array that needs to contains the following fields: 
%   radius (Rad)        Radii of the cylinders, vector
%   length (Len)        Lengths of the cylinders, vector
%   parent (CPar)       Parents of the cylinders, vector
% inputs.GrowthVolFac   The factor "fac", defines the upper and lower
%                         allowed radius from the predicted one:
%                         1/fac*predicted_rad <= rad <= fac*predicted_rad
% ---------------------------------------------------------------------

% Changes from version 1.0.0 to 2.0.0, 16 Sep 2021:
% 1) Changed the roles of RADIUS and GROWTH_VOLUME in the allometry, i.e.
%    the radius is now predicted from the growth volume
% 2) Do not increase the radius of the branch tip cylinders 

disp('----------')
disp('Growth volume based correction of cylinder radii:')

Rad = double(cylinder.radius);
Rad0 = Rad;
Len = double(cylinder.length);
CPar = cylinder.parent;
CExt = cylinder.extension;

initial_volume = round(1000*pi*sum(Rad.^2.*Len));
disp([' Initial_volume (L): ',num2str(initial_volume)])

%% Define the child cylinders for each cylinder
n = length(Rad);
CChi = cell(n,1);
ind = (1:1:n)';
for i = 1:n
  CChi{i} = ind(CPar == i);
end

%% Compute the growth volume
GrowthVol = zeros(n,1); % growth volume
S = cellfun('length',CChi);
modify = S == 0;
GrowthVol(modify) = pi*Rad(modify).^2.*Len(modify);
parents = unique(CPar(modify));
if parents(1) == 0
  parents = parents(2:end);
end
while ~isempty(parents)
  V = pi*Rad(parents).^2.*Len(parents);
  m = length(parents);
  for i = 1:m
    GrowthVol(parents(i)) = V(i)+sum(GrowthVol(CChi{parents(i)}));
  end
  parents = unique(CPar(parents));
  if parents(1) == 0
    parents = parents(2:end);
  end
end

%% Fit the allometry: Rad = a*GV^b;
options = optimset('Display','off');
X = lsqcurvefit(@allometry,[0.5 0.5 0],GrowthVol,Rad,[],[],options);
disp(' Allometry model parameters R = a*GV^b+c:')
disp(['   Multiplier a: ', num2str(X(1))])
disp(['   Exponent b: ', num2str(X(2))])
if length(X) > 2
disp(['   Intersect c: ', num2str(X(3))])
end

%% Compute the predicted radius from the allometry
PredRad = allometry(X,GrowthVol);

%% Correct the radii based on the predictions
% If cylinder's radius is over fac-times or under 1/fac-times the
% predicted radius, then correct the radius to match the allometry
fac = inputs.GrowthVolFac;
modify = Rad < PredRad/fac | Rad > fac*PredRad;
modify(Rad < PredRad/fac & CExt == 0) = 0; % Do not increase the radius at tips
CorRad = PredRad(modify);

% Plot allometry and radii modification
gvm = max(GrowthVol);
gv = (0:0.001:gvm);
PRad = allometry(X,gv);
figure(1)
plot(GrowthVol,Rad,'.b','Markersize',2)
hold on
plot(gv,PRad,'-r','Linewidth',2)
plot(gv,PRad/fac,'-g','Linewidth',2)
plot(gv,fac*PRad,'-g','Linewidth',2)
hold off
grid on
xlabel('Growth volume (m^3)')
ylabel('Radius (m)')
legend('radius','predicted radius','minimum radius','maximum radius','Location','NorthWest')

figure(2)
histogram(CorRad-Rad(modify))
xlabel('Change in radius')
title('Number of cylinders per change in radius class')

% Determine the maximum radius change
R = Rad(modify);
D = max(abs(R-CorRad)); % Maximum radius change
J = abs(R-CorRad) == D;
D = CorRad(J)-R(J);

% modify the radius according to allometry
Rad(modify) = CorRad; 
cylinder.radius = Rad;

disp([' Modified ',num2str(nnz(modify)),' of the ',num2str(n),' cylinders'])
disp([' Largest radius change (cm): ',num2str(round(1000*D)/10)])
corrected_volume = round(1000*pi*sum(Rad.^2.*Len));
disp([' Corrected volume (L): ', num2str(corrected_volume)])
disp([' Change in volume (L): ', num2str(corrected_volume-initial_volume)])
disp('----------')

% % Plot cylinder models where the color indicates change (green = no change, 
% % red = decreased radius, cyan = increased radius)
% cylinder.branch = ones(n,1);
% cylinder.BranchOrder = ones(n,1);
% I = Rad < Rad0;
% cylinder.BranchOrder(I) = 2;
% I = Rad > Rad0;
% cylinder.BranchOrder(I) = 3;
% plot_cylinder_model(cylinder,'order',3,20,1)
% 
% cyl = cylinder;
% cyl.radius = Rad0;
% plot_cylinder_model(cyl,'order',4,20,1)

end % End of main function


function F = allometry(x,xdata)
F = x(1)*xdata.^x(2)+x(3);
end
