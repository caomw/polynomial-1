% Author: Mathieu Salzmann, NICTA.
%
function Mesh = LoadMesh(PtsName, TriName)

Mesh.coords = load(PtsName);
nPts = size(Mesh.coords,1);
Mesh.neighbors = zeros(nPts,nPts);
trig = load(TriName);

for i=1:size(trig,1)
    Mesh.neighbors(trig(i,1)+1,trig(i,2)+1) = 1;
    Mesh.neighbors(trig(i,2)+1,trig(i,1)+1) = 1;
    Mesh.neighbors(trig(i,1)+1,trig(i,3)+1) = 1;
    Mesh.neighbors(trig(i,3)+1,trig(i,1)+1) = 1;
    Mesh.neighbors(trig(i,2)+1,trig(i,3)+1) = 1;
    Mesh.neighbors(trig(i,3)+1,trig(i,2)+1) = 1;
end

Mesh.E = ComputeEdges(Mesh);
