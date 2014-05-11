% Author: Mathieu Salzmann, NICTA.
%
% Exhaustively evaluate all possible combinations of the slaves to find the
% best ones. This is only done when there are a small enough number of
% discrete solutions for each slave.
function [x,f,minerr] = find_best_combination(x,f,xall,fall,indsgx,ng,nxpg,dimx,nx,bname,groups,vars)

minerr = inf;
if(length(xall)==1)
    for i=1:length(xall{1})
        xval = zeros(ng*nxpg,dimx);
        for j=1:ng-1
            xval((j-1)*nxpg+1:j*nxpg,:) = x{j};
        end
        xval((ng-1)*nxpg+1:ng*nxpg,:) = xall{1}{i};
        xt = zeros(nx,dimx);
        for j=1:nx
            xt(j,:) = mean(xval(indsgx{j},:));
        end
        err = feval(bname,xt,groups,vars);
        if(err<minerr)
            minerr = err;
            indbest = i;
        end
    end
    x{length(x)+1} = xall{1}{indbest};
    f{length(f)+1} = fall{1}{indbest};
else
    for i=1:length(xall{1})
        tmpx = x;
        tmpx{length(tmpx)+1} = xall{1}{i};
        tmpf = f;
        tmpf{length(tmpf)+1} = fall{1}{i};
        [xtmp,ftmp,err] = find_best_combination(tmpx,tmpf,xall(2:end),fall(2:end),indsgx,ng,nxpg,dimx,nx,bname,groups,vars);
        if(err<minerr)
            minerr = err;
            xbest = xtmp;
            fbest = ftmp;
        end
    end
    x = xbest;
    f = fbest;
end

