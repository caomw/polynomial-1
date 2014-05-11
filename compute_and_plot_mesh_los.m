% Author: Mathieu Salzmann, NICTA.
%
% Compute the 3D mesh from the depths x and the lines of sight Q.
% Plot the resulting mesh.
function compute_and_plot_mesh_los(x,vars,Mesh)

[Q,E,Eg,W] = deal(vars{:});

y = repmat(x',3,1).*Q;

figure(1);
clf;
PlotMesh(Mesh,'b',1);
Mesh.coords = y';
PlotMesh(Mesh,'r',1);
view(-25,25);
pause(0.1);

