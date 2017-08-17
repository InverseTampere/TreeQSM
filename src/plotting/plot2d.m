function h = plot2d(X,Y,fig,strtit,strx,stry,leg,E)

% 2D-plots, where the data (X and Y), figure number, title, xlabel, ylabel, 
% legends and error bars can be specied with the inputs.

lw = 1.5; % linewidth
n = size(Y,1);
if n < 9
    col = ['-b '; '-r '; '-g '; '-c '; '-m '; '-k '; '-y '; '-.b'];
else
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
    if n > 20
        k = ceil(n/20);
        col = repmat(col,[k 1]);
    end
end
figure(fig)
if nargin <= 7
    % plots without errorbars
    if ~iscell(Y)
        if ~isempty(X)
            if n < 9
                h = plot(X(1,:),Y(1,:),'-b','Linewidth',lw);
            else
                h = plot(X(1,:),Y(1,:),'Color',col(1,:),'Linewidth',lw);
            end
        else
            if n < 9
                h = plot(Y(1,:),'-b','Linewidth',lw);
            else
                h = plot(Y(1,:),'Color',col(1,:),'Linewidth',lw);
            end
        end
        if n > 1
            hold on
            if ~isempty(X)
                if n < 9
                    for i = 2:n
                        plot(X(i,:),Y(i,:),col(i,:),'Linewidth',lw)
                    end
                else
                    for i = 2:n
                        plot(X(i,:),Y(i,:),'Color',col(i,:),'Linewidth',lw)
                    end
                end
            else
                if n < 9
                    for i = 2:n
                        plot(Y(i,:),col(i,:),'Linewidth',lw)
                    end
                else
                    for i = 2:n
                        plot(Y(i,:),'Color',col(i,:),'Linewidth',lw)
                    end
                end
            end
            hold off
        end
    else
        if ~isempty(X)
            x = X{1};
        end
        y = Y{1};
        if ~isempty(X)
            if n < 9
                h = plot(x,y,'-b','Linewidth',lw);
            else
                h = plot(x,y,'Color',col(1,:),'Linewidth',lw);
            end
        else
            if n < 9
                h = plot(y,'-b','Linewidth',lw);
            else
                h = plot(y,'Color',col(1,:),'Linewidth',lw);
            end
        end
        if n > 1
            hold on
            if ~isempty(X)
                for i = 2:n
                    x = X{i};
                    y = Y{i};
                    if n < 9
                        plot(x,y,col(i,:),'Linewidth',lw)
                    else
                        plot(x,y,'Color',col(i,:),'Linewidth',lw)
                    end
                end
            else
                for i = 2:n
                    y = Y{i};
                    if n < 9
                        plot(y,col(i,:),'Linewidth',lw)
                    else
                        plot(y,'Color',col(i,:),'Linewidth',lw)
                    end
                end
            end
            hold off
        end
    end
    
else
    % plots with errorbars
    if ~iscell(Y)
        if ~isempty(X)
            if n < 9
                h = errorbar(X(1,:),Y(1,:),E(1,:),'-b','Linewidth',lw);
            else
                h = errorbar(X(1,:),Y(1,:),E(1,:),'Color',col(1,:),'Linewidth',lw);
            end
        else
            if n < 9
                h = errorbar(Y(1,:),E(1,:),'-b','Linewidth',lw);
            else
                h = errorbar(Y(1,:),E(1,:),'Color',col(1,:),'Linewidth',lw);
            end
        end
        if n > 1
            hold on
            if ~isempty(X)
                if n < 9
                    for i = 2:n
                        errorbar(X(i,:),Y(i,:),E(1,:),col(i,:),'Linewidth',lw)
                    end
                else
                    for i = 2:n
                        errorbar(X(i,:),Y(i,:),E(1,:),'Color',col(i,:),'Linewidth',lw)
                    end
                end
            else
                if n < 9
                    for i = 2:n
                        errorbar(Y(i,:),E(1,:),col(i,:),'Linewidth',lw)
                    end
                else
                    for i = 2:n
                        errorbar(Y(i,:),E(1,:),'Color',col(i,:),'Linewidth',lw)
                    end
                end
            end
            hold off
        end
    else
        if ~isempty(X)
            x = X{1};
        end
        y = Y{1};
        e = E{1};
        if ~isempty(X)
            if n < 9
                h = errorbar(x,y,e(1,:),'-b','Linewidth',lw);
            else
                h = errorbar(x,y,'Color',e(1,:),col(1,:),'Linewidth',lw);
            end
        else
            if n < 9
                h = errorbar(y,e(1,:),'-b','Linewidth',lw);
            else
                h = errorbar(y,e(1,:),'Color',col(1,:),'Linewidth',lw);
            end
        end
        if n > 1
            hold on
            if ~isempty(X)
                for i = 2:n
                    x = X{i};
                    y = Y{i};
                    e = E{i};
                    if n < 9
                        h = errorbar(x,y,e(1,:),col(i,:),'Linewidth',lw);
                    else
                        h = errorbar(x,y,e(1,:),'Color',col(i,:),'Linewidth',lw);
                    end
                end
            else
                for i = 2:n
                    y = Y{i};
                    e = E{i};
                    if n < 9
                        errorbar(y,e(1,:),col(i,:),'Linewidth',lw);
                    else
                        errorbar(y,e(1,:),'Color',col(i,:),'Linewidth',lw);
                    end
                end
            end
            hold off
        end
    end
end

grid on
t = title(strtit);
x = xlabel(strx);
y = ylabel(stry);
if nargin > 6
    legend(leg)
end
set(gca,'fontsize',12)
set(gca,'FontWeight','bold')
set(t,'fontsize',12)
set(t,'FontWeight','bold')
set(x,'fontsize',12)
set(x,'FontWeight','bold')
set(y,'fontsize',12)
set(y,'FontWeight','bold')