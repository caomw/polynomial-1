% Author: Mathieu Salzmann, NICTA.
%
% Main algorithms: DD-Poly (if rho1 = 0) or ADMM-Poly
% 
% INPUT:
% groups: indices of the different slaves
% dimx: dimension of each individual variables
% fname: function that solves a slave problem
% niter: maximum number of iterations of dual decomposition
% grad_rule: subgradient update rule: (1) adaptive, (0) non-summable diminishing step length
% from eta0 to gamma: parameters, see in the paper
% vars: additional values for the problem (e.g., lines of sight)
% bname: function that compute the global energy
% varargin: values to potentially plot the result at each iteration
%
% OUTPUT:
% xbest: final global solution (either by convergence, or max number of iterations)
% xallbest: corresponding solutions for each slave
% primal: primal energy at each iteration
% dual: dual energy at each iteration
function [xbest,xallbest,primal,dual] = dd_admm_poly(groups,dimx,fname,niter,grad_rule,eta0,delta0,tau0,tau1,t0,rho0,rho1,gamma,vars,bname,varargin)

% Different global and group variables
[ng,nxpg] = size(groups);
nx = max(max(groups));
xbest = zeros(nx,dimx);
xm = zeros(nx,dimx);
xall = zeros(ng*nxpg,dimx);
xallbest = xall;
xtmp = cell(ng,1);

% Cells that will contain the funtion values and gradients for the different slaves
ftmp = cell(ng,1);
gpart = cell(ng,1);
gradstr = cell(ng,1);
for i=1:ng
    gpart{i} = sym([]);
    gradstr{i} = [];
end
nsol = zeros(ng,1);

% Cell for group indexing
indsgx = cell(nx,1);
for i=1:nx
    indsgx{i} = [];
end
for i=1:ng
    for j=1:nxpg
        indsgx{groups(i,j)} = [indsgx{groups(i,j)};(i-1)*nxpg+j];
    end
end

% Lagrange multipliers
lambda = zeros(size(xall));

rgap = 1; % relative primal-dual gap
agap = 1; % absolute primal-dual gap

bestdual = -inf;
bestprimal = inf;
primal = zeros(1,niter);
dual = zeros(1,niter);

it = 1; % iteration number
delta = delta0;

while((it<=niter)&&(rgap>1e-3)&&(agap>1e-6))
    fprintf('Iteration %d\n', it);
    % Compute the augmented Lagrangian weight
    if(it==1)
        rho = rho0;
    else
        rho = rho0 + (rho1-rho0)/(1 + exp(-gamma*(it-t0)));
    end
    
    % Solve the slaves
    parfor i=1:ng
        % one folder per slave
        dirname = sprintf('%d',i);
        cd(dirname)
        [xtmp{i},ftmp{i},gpart{i},gradstr{i},nsol(i)] = feval(fname,i,groups(i,:),lambda((i-1)*nxpg+1:i*nxpg,:),gpart{i},gradstr{i},rho,xm,vars);
        cd ..
    end
    
    if(isempty(find(nsol==0, 1))) % Everything went well, we have at least one solution for each slave
        
        ntot = sum(nsol);
        if((ntot>ng)&& (ntot<1024)) % Some groups have more than one solution, but there is a reasonable total number of solution
            [xtmp,ftmp] = find_best_combination([],[],xtmp,ftmp,indsgx,ng,nxpg,dimx,nx,bname,groups,vars); % Let's try and find the best combination
        else % There is a single solution per slave, or too many solutions
            % Let's take the first solution for each slave
            for i=1:ng
                xtmp{i} = xtmp{i}{1};
                ftmp{i} = ftmp{i}{1};
            end
        end

        % Compute the slave and global solutions, as well as primal and
        % dual function values.
        for i=1:ng
            xall((i-1)*nxpg+1:i*nxpg,:) = xtmp{i};
            dual(it) = dual(it) + ftmp{i};
        end
        for i=1:nx
            xm(i,:) = mean(xall(indsgx{i},:),1);
        end
        primal(it) = feval(bname,xm,groups,vars);
        if(primal(it)<=bestprimal)
            bestprimal = primal(it);
            xallbest = xall;
        end
        improve = dual(it)-bestdual;
        if(dual(it)>=bestdual)
            bestdual = dual(it);
        end
                
        % Compute the projected subgradient
        gl = xall;
        for i=1:nx
           gl(indsgx{i},:) = gl(indsgx{i},:) - repmat(xm(i,:),length(indsgx{i}),1);
        end
        
        % Compute the subgradient step
        ndg = norm(reshape(gl,[],1));
        approx = bestdual + delta;
        if(improve>=delta)
            delta = tau0*delta;
        else
            delta = max(tau1*delta,1e-3*delta0);
        end
        if(grad_rule>0)
            eta = (approx-dual(it))/(ndg^2);
        else
            eta = (eta0/(1 + it/t0))/ndg;
        end
        
        % Update the Lagrange multipliers
        lambda = lambda + eta*gl;
        
        % Compute the gaps
        agap = abs(primal(it)-dual(it));
        rgap = agap/primal(it);

        % Print some statistics
        fprintf('It %d: primal = %d, dual = %d, rel gap = %d\n',it,primal(it),dual(it),rgap);
        
        % Plot the results if needed
        if((length(varargin)>=2) && varargin{1})
            feval(varargin{2},xm,vars,varargin{3:end});
        end
        
    else % Something bad happened. Let's just take a smaller subgradient step and start over.
        fprintf('No solution found for at least one slave!\n');
        lambda = lambda - eta*gl;
        eta = 0.5*eta;
        lambda = lambda + eta*gl;
        dual(it) = dual(it-1);
        primal(it) = primal(it-1);
    end
    it = it+1;
end

primal = primal(1:it-1);
dual = dual(1:it-1);
for i=1:nx
    xbest(i,:) = mean(xallbest(indsgx{i},:),1);
end
