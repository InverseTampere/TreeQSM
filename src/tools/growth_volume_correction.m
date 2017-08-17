function cylinder = growth_volume_correction(cylinder,inputs)

% Use Jan Hackenberg's growth volume allometry approach to modify the radius
% of too large and too small cylinders. Uses the allometry 
%       GrowthVolume = a * Radius^b + c
% If cylinder's growth volume is over fac-times or under 1/fac-times the 
% predicted growth volume, then correct the radius to such that it 
% corresponds to the allometry.
% More details can be from Jan's "SimpleTree" papers and documents.

Rad = double(cylinder.radius);
Len = double(cylinder.length);
CExt = cylinder.extension;
CChi = cylinder.ChildCyls;
CPar = cylinder.parent;

%initial_volume = round(1000*pi*sum(Rad.^2.*Len));
%disp(initial_volume)

% Modify CChi, child cylinders, to include the extension for easier growth
% volume computation
n = length(Rad);
for i = 1:n
    if CExt(i) > 0
        CChi{i} = [CChi{i}; CExt(i)];
    end
end

% Compute the growth volume
GrowthVol = zeros(n,1); % growth volume
S = cellfun('length',CChi);
I = S == 0;
GrowthVol(I) = pi*Rad(I).^2.*Len(I);
parents = unique(CPar(I));
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

% Fit the allometry: GV = a*Rad^b+c;
options = optimset('Display','off');
X = lsqcurvefit(@allometry,[100 2.5 0],Rad,GrowthVol,[],[],options);
%disp(X(1))
%disp(X(2))
%disp(X(3))

% Compute the predicted growth volume from the allometry
PredGrowthVol = allometry(X,Rad);

% If cylinder's growth volume is over fac-times or under 1/fac-times the 
% predicted growth volume, then correct the radius to such that it 
% corresponds to the allometry
fac = inputs.GrowthVolFac;
I = GrowthVol < PredGrowthVol/fac | GrowthVol > fac*PredGrowthVol;
R = exp(log((GrowthVol(I)-X(3))/X(1))/X(2));
R = real(R);

% % Plot allometry and radii modification 
% rm = max(Rad);
% r = (0:0.001:rm);
% PGV = allometry(X,r);
% figure(1)
% plot(Rad,GrowthVol,'.b')
% hold on
% plot(r,PGV,'-r')
% plot(r,PGV/fac,'-g')
% plot(r,fac*PGV,'-g')
% hold off
% grid on
% 
% a = max(max(Rad(I)),max(R));
% figure(3)
% plot(Rad(I),R,'.b')
% hold on
% plot(Rad(~I),Rad(~I),'og')
% plot([0; a],[0; a],'-r')
% hold off
% axis equal
% grid on
%
% figure(4)
% histogram(Rad(I)-R,100)
% D = max(abs(Rad(I)-R));

Rad(I) = R; % modify the radius according to allometry
cylinder.radius = Rad;

%disp([n nnz(I)])
%disp(D)
%corrected_volume = round(1000*pi*sum(Rad.^2.*Len));
%disp(corrected_volume)
%disp(initial_volume-corrected_volume)
%pause

end


function F = allometry(x,xdata)
    F = x(1)*xdata.^x(2)+x(3);    
end