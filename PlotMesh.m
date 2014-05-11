% Author: Mathieu Salzmann, NICTA.
%
function PlotMesh(Mesh,mycolor,figid)

if(nargin<3)
    figid = 1;
end
figure(figid);

hold on;

if nargin<2
    mycolor = 'k';
end

daspect('manual');
daspect([1 1 1]);
pbaspect('manual');
pbaspect([1 1 1]);

for i=1:size(Mesh.coords,1)
    for j=i+1:size(Mesh.coords,1)
        if Mesh.neighbors(i,j)==1
            plot3([Mesh.coords(i,1),Mesh.coords(j,1)], [Mesh.coords(i,2),Mesh.coords(j,2)], [Mesh.coords(i,3),Mesh.coords(j,3)],mycolor,'linewidth',2);
        end
    end
end

hold off;
