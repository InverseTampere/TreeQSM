function plot_cone_model(cylinder,fig,nf,alp,Ind)

% Plots the given cylinder model as truncated cones defined by the cylinders.
% cylinder  Structure array containin the cylinder info 
%               (radius, length, start, axis, BranchOrder)
% fig       Figure number
% nf        Number of facets in the cyliders (in the thickest cylinder, 
%               scales down with radius to 4 which is the minimum)
% alp       Alpha value (1 = no trasparency, 0 = complete transparency)
% Ind       Indexes of cylinders to be plotted from a subset of cylinders
%               (Optional, if not given then all cylinders are plotted)


if isstruct(cylinder)
    Rad = cylinder.radius;
    Len = cylinder.length;
    Sta = cylinder.start;
    %Sta = mat_vec_subtraction(Sta,Sta(1,:));
    Axe = cylinder.axis;
    Bran = cylinder.branch;
    PiB = cylinder.PositionInBranch;
    nb = max(Bran);
else
    Rad = cylinder(:,1);
    Len = cylinder(:,2);
    Sta = cylinder(:,3:5);
    %Sta = mat_vec_subtraction(Sta,Sta(1,:));
    Axe = cylinder(:,6:8);
    Bran = cylinder(:,12);
    PiB = cylinder(:,14);
    nb = max(Bran);
end
if nargin == 5
    Rad = Rad(Ind);
    Len = Len(Ind);
    Sta = Sta(Ind,:);
    Axe = Axe(Ind,:);
end

nc = size(Rad,1);

Cir = cell(nf,2);
for i = 4:nf
    Cir{i,1} = [cos((1/i:1/i:1)*2*pi)' sin((1/i:1/i:1)*2*pi)' zeros(i,1)];
    Cir{i,2} = [(1:1:i)' (i+1:1:2*i)' [(i+2:1:2*i)'; i+1] [(2:1:i)'; 1]];
end

Vert = zeros(2*nc*(nf+1),3);
Facets = zeros(nc*(nf+1),4);
t = 1;
f = 1;

% Scale, rotate and translate the standard cylinders
Ind = (1:1:nc)';
for j = 1:nb
    I = Bran == j;
    I = Ind(I);
    if ~isempty(I)
        P = PiB(I);
        [P,J] = sort(P);
        I = I(J);
        n = ceil(sqrt(mean(Rad(I))/Rad(1))*nf);
        n = min(n,nf);
        n = max(n,4);
        C0 = Cir{n,1};
        m = length(I);
        for i = 1:m
            C = C0;
            
            % Scale radius
            C(1:n,1:2) = Rad(I(i))*C(1:n,1:2);
            if i == m
                % Define the last circle of the branch
                C1 = C;
                C1(:,1:2) = min(0.005/Rad(I(i)),1)*C(:,1:2);
            end
            
            % Rotate
            if i == 1
                ang = real(acos(Axe(I(i),3)));
                Axis = cross([0 0 1]',Axe(I(i),:)');
                Rot = rotation_matrix(Axis,ang);
                C = C*Rot';
            elseif i > 1
                ang = real(acos(Axe(I(i),3)));
                Axis = cross([0 0 1]',Axe(I(i),:)');
                Rot = rotation_matrix(Axis,ang);
                C = C*Rot';
                %%% Should be somehow corrected so that high angles between
                %%% cylinders do not cause narrowing the surface!!!
                
                
                if i == m
                    ang = real(acos(Axe(I(i),3)));
                    Axis = cross([0 0 1]',Axe(I(i),:)');
                    Rot = rotation_matrix(Axis,ang);
                    C1 = C1*Rot';
                end
            end
            
            % Translate
            C = mat_vec_subtraction(C,-Sta(I(i),:));
            if i == m
                C1 = mat_vec_subtraction(C1,-(Sta(I(i),:)+Len(I(i))*Axe(I(i),:)));
            end
            
            % Save the new vertices
            Vert(t:t+n-1,:) = C;
            if i == m
                t = t+n;
                Vert(t:t+n-1,:) = C1;
            end
            t = t+n;
            
            % Define the new facets
            if i == 1 && i == m
                Facets(f:f+n-1,:) = Cir{n,2}+t-2*n-1;
                f = f+n;
            elseif i > 1 && i < m 
                Facets(f:f+n-1,:) = Cir{n,2}+t-2*n-1;
                f = f+n;
            elseif i > 1 && i == m
                Facets(f:f+n-1,:) = Cir{n,2}+t-3*n-1;
                f = f+n;
                Facets(f:f+n-1,:) = Cir{n,2}+t-2*n-1;
                f = f+n;
            end
        end
    end
end

t = t-1;
f = f-1;
Vert = Vert(1:t,:);
Facets = Facets(1:f,:);
fvd = [139/255*ones(f,1) 69/255*ones(f,1) 19/255*ones(f,1)];

figure(fig)
plot3(Vert(1,1),Vert(1,2),Vert(1,3))
patch('Vertices',Vert,'Faces',Facets,'FaceVertexCData',fvd,'FaceColor','flat')
alpha(alp)
axis equal
grid on
view(-37.5,30)
