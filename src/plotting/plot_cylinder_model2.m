function plot_cylinder_model2(cylinder,fig,nf,alp,Ind)

% Plots the cylinder model.
% cylinder  Structure array containin the cylinder info 
%               (radius, length, start, axis, BranchOrder)
% fig       Figure number
% nf        Number of facets in the cyliders (in the thickest cylinder, 
%               scales down with radius to 4 which is the minimum)
% alp       Alpha value (1 = no trasparency, 0 = complete transparency)
% Ind       Indexes of cylinders to be plotted from a subset of cylinders
%               (Optional, if not given then all cylinders are plotted)


Rad = cylinder.radius;
Rad2 = cylinder.TopRadius;
Len = cylinder.length;
Sta = cylinder.start;
Sta = mat_vec_subtraction(Sta,Sta(1,:));
Axe = cylinder.axis;
BOrd = cylinder.BranchOrder;
if nargin == 5
    Rad = Rad(Ind);
    Len = Len(Ind);
    Sta = Sta(Ind,:);
    Axe = Axe(Ind,:);
    BOrd = BOrd(Ind);
end

nc = size(Rad,1);

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

N = max(BOrd)+1;
if N <= 20
    col = col(1:N,:);
else
    m = ceil(N/20);
    col = repmat(col,[m,1]);
    col = col(1:N,:);
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
    m = size(C,1);
    C(1:m/2,1:2) = Rad(i)*C(1:m/2,1:2);
    C(m/2+1:m,1:2) = Rad2(i)*C(m/2+1:m,1:2);
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
    fvd(f:f+n-1,:) = repmat(col(BOrd(i)+1,:),[n 1]);
    t = t+2*n;
    f = f+n;
end
t = t-1;
f = f-1;
Vert = Vert(1:t,:);
Facets = Facets(1:f,:);
fvd = fvd(1:f,:);

figure(fig)
plot3(Vert(1),Vert(2),Vert(3))
patch('Vertices',Vert,'Faces',Facets,'FaceVertexCData',fvd,'FaceColor','flat')
alpha(alp)
axis equal
grid on
view(-37.5,30)
