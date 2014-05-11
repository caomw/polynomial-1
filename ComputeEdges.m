% Author: Mathieu Salzmann, NICTA.
%
function E = ComputeEdges(Mesh)

pts = Mesh.coords;
Adj = Mesh.neighbors;
nPts = size(pts,1);

E = zeros(0,3);
ind = 1;
for i=1:nPts
    for j=i+1:nPts
        if(Adj(i,j)==1)
            E(ind,1:2) = [i,j];
            t = pts(j,:) - pts(i,:);
            E(ind,3) = sqrt(t*t');
            ind = ind+1;
        end
    end
end

