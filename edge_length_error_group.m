% Author: Mathieu Salzmann, NICTA.
%
% Compute the solution to one specific slave.
% Compute the gradient of the slave energy, write the corresponding HOM4PS2
% file, solve the system of polynomial equations with HOM4PS2 and read the
% solutions and find the globally optimal one.
%
% INPUT:
% subInd: slave index
% vinds: indices of the vertices involved in the slave
% lambda: Lagrange multipliers
% gpart: Energy function for the slave (computed during the first call)
% gradstr: string representation of the gradient of the energy
% rho: weight of the augmented Lagrangian
% xp: current global solution in ADMM
% vars: additional variables (lines of sight, edges, slave edges, weights)
%
% OUTPUT:
% d: depths of the mesh vertices
% fv: energy value
% gpart and gradstr: see above
% nsol: number of discrete solutions (ideally 1, but disambiguated later if
% several)
function [d,fv,gpart,gradstr,nsol] = edge_length_error_group(subInd,vinds,lambda,gpart,gradstr,rho,xp,vars)

[Q,E,Eg,W] = deal(vars{:});
nvar = length(vinds);

% Compute the energy and its gradient and save it in string form to make it easier to
% interact with HOM4PS2.
if(isempty(gpart))
    for i=1:nvar
       name = sprintf('x%d',vinds(i));
       x(vinds(i)) = sym(name,'real');
    end
    x = x';
    for i=1:length(Eg(subInd,:))
        tmp = x(E(Eg(subInd,i),1))*Q(:,E(Eg(subInd,i),1)) - x(E(Eg(subInd,i),2))*Q(:,E(Eg(subInd,i),2));
        f(i) = tmp'*tmp - E(Eg(subInd,i),3)^2;
    end
    gpart = 0.5*((W(subInd,:).*f)*(f'))/size(E,1);
    gradstr = cell(nvar,1);
    for i=1:nvar
        gradstr{i} = [];
        name = sprintf('x%d',vinds(i));
        gradpart = diff(gpart,name);
        [cf,tf] = coeffs(gradpart,x(vinds));
        tmp = char(vpa(cf(1)));
        if(tmp(1)=='-')
            gradstr{i} = [gradstr{i},['- (' tmp(2:end) ')*']];
        else
            gradstr{i} = [gradstr{i},['(' tmp ')*']];
        end
        gradstr{i} = [gradstr{i},['(' char(tf(1)) ')']];
        for j=2:length(cf)
            tmp = char(vpa(cf(j)));
            if(tmp(1)=='-')
                gradstr{i} = [gradstr{i},[' - (' tmp(2:end) ')*']];
            else
                gradstr{i} = [gradstr{i},[' + (' tmp ')*']];
            end
            gradstr{i} = [gradstr{i},['(' char(tf(j)) ')']];
        end
    end
end

% Write the file that will act as input to HOM4PS2.
name = sprintf('mesh3x3_%d',subInd);
fd = fopen(name,'w');
fprintf(fd,'{\n');
for i=1:nvar
    fprintf(fd,'%s',gradstr{i});
    if(rho>1e-3)
        fprintf(fd,' + (%f)*(x%d)',rho,vinds(i));
        if(xp(vinds(i))>0)
            fprintf(fd,' - %f',rho*xp(vinds(i)));
        elseif(xp(vinds(i))<0)
            fprintf(fd,' + %f',-rho*xp(vinds(i)));
        end
    end
    if(lambda(i)>0)
        fprintf(fd,' + %f;\n',lambda(i));
    elseif(lambda(i)<0)
        fprintf(fd,' %f;\n',lambda(i));
    else
        fprintf(fd,';\n');
    end
end
fprintf(fd,'}\n');
fclose(fd);

% Find the solutions of the polynomial system with HOM4PS2.
fprintf('Calling hom4ps2...\n');
tmp = pwd;
tmpn = sprintf('./hom4ps2 %s/%s 1>NUL',tmp,name);
system(tmpn);

% Read all the solutions obtained by HOM4PS2 and find the best one (or best ones
% if there is a tie between multiple solutions).
d{1} = zeros(nvar,1);
fv{1} = inf;
nsol = 0;
sols = [];
fd = fopen('data.roots','r');
if(fd)
    c = fscanf(fd,'%c',1);
    while(c == '(')
        tmp = zeros(1,nvar);
        test = 1;
        for j=1:nvar
            tmp(j) = fscanf(fd,'%f',1);
            if(tmp(j)<1) % no negative depth
                test = 0;
            end
            fscanf(fd,'%c',2);
            im = fscanf(fd,'%f',1);
            if(abs(im) > 1e-6) % no imaginary solutions
                test = 0;
            end
            fscanf(fd,'%c',3);
        end
        if(test)
            sols = [sols;tmp];
        end
        fgetl(fd);
        fgetl(fd);
        fgetl(fd);
        fgetl(fd);
        c = fscanf(fd,'%c',1);
    end
    fgetl(fd);
    indvars = zeros(1,nvar);
    fscanf(fd,'%c',2);
    for j=1:nvar
        tmp = fscanf(fd,'%d',1);
        indvars(j) = find(vinds==tmp);
        fscanf(fd,'%c',3);
    end
    fclose(fd);
    if(~isempty(sols))
        sols(:,indvars) = sols;
        for i=1:size(sols,1)
            xtmp(vinds) = sols(i,:);
            xtmp = xtmp';
            fp = zeros(1,length(Eg(subInd,:)));
            for j=1:length(Eg(subInd,:))
                tmp = xtmp(E(Eg(subInd,j),1))*Q(:,E(Eg(subInd,j),1)) - xtmp(E(Eg(subInd,j),2))*Q(:,E(Eg(subInd,j),2));
                fp(j) = tmp'*tmp - E(Eg(subInd,j),3)^2;
            end
            fpart = 0.5*((W(subInd,:).*fp)*(fp'))/size(E,1);
            ftmp = fpart + lambda'*sols(i,:)';
            if(abs(ftmp-fv{1})<1e-6) % we have a tie
                nsol = nsol+1;
                fv{nsol} = ftmp;
                d{nsol} = sols(i,:)';
            elseif(ftmp<fv{1})
                fv = [];
                fv{1} = ftmp;
                fov = [];
                fov{1} = fpart;
                d = [];
                d{1} = sols(i,:)';
                nsol = 1;
            end
        end
    end
else
    fprintf('Error: no solution was found.\n');
end

