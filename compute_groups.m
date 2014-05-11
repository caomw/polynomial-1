% Author: Mathieu Salzmann, NICTA.
%
% Compute the groups for dual decomposition.
% A wxh mesh is decomposed into small 2x2 overlapping square meshes.
% G: matrix of original indices for each group
% Eg: Mesh edges for each group
% W: Weights to account for the fact that some edges appear in multiple
% groups (this enforces having the same energy as before despite the
% decomposition).
function [G,Eg,W] = compute_groups(w,h,E)

G = zeros((w-1)*(h-1),4);
Eg = zeros((w-1)*(h-1),5);
W = zeros((w-1)*(h-1),5);

for i=1:h-1
    for j=1:w-1
        G((i-1)*(w-1)+j,:) = [(i-1)*w+j,(i-1)*w+j+1,i*w+j,i*w+j+1];
        ind = 1;
        for k=1:size(E,1)
            i1 = find(E(k,1)==G((i-1)*(w-1)+j,:), 1);
            i2 = find(E(k,2)==G((i-1)*(w-1)+j,:), 1);
            if(~isempty(i1) && ~isempty(i2))
                Eg((i-1)*(w-1)+j,ind) = k;
                ind = ind + 1;
            end
        end
    end
end

for k=1:size(E,1)
    inds = find(Eg(:)==k);
    W(inds) = 1/length(inds);
end
