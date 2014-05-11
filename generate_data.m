% Author: Mathieu Salzmann, NICTA.
%
% This code performs non-rigid 3D reconstruction of a 3x3 square mesh, as
% discussed in Section 4.2 of my CVPR 2013 paper "Continuous Inference in
% Graphical Models with Polynomial Energies". Reconstruction is performed 
% for 100 random deformations.

%% Add the folder to the path
addpath(pwd); % Assumes that your current folder is the one where the code is

%% Parameter settings

NDefs = 100; % Number of deformations to reconstruct (up to 100)
noise = 2; % Std of the image location (input) noise

finaldisplay = 1; % Display the final result
display = 0; % Display the result at each iteration of dual decomposition
niter = 100; % Maximum number of iterations of dual decomposition
grad_rule = 0; % Adaptive rule (1) or non-summable diminishing step length (0)
eta0 = 1e1; % See Eq. 8
delta0 = eta0; % See Eq. 7
rho1 = 1e2; % Final augmented Lagrangian weight in ADMM-Poly. If 0, then DD-Poly is used.

%% Load data and generate noisy image locations

basename = 'data/';

Ygt = load([basename 'mesh_3x3_100_deformations_ground_truth.txt']); % Ground-truth 3D meshes

Ref = LoadMesh([basename 'mesh_3x3.pts'],[basename 'mesh_3x3.tri']); % Load reference mesh
Mesh = Ref;
nPts = size(Ref.coords,1);
E = Ref.E; % Mesh edges

int = load([basename 'cam.intr']); % Camera intrinsic parameters
ext = load([basename 'cam.ext']); % Camera extrinsic parameters

% Project the 3D points of each mesh to image locations
uv = zeros(100,2*nPts);
for i=1:NDefs
    tmp = int*ext*[reshape(Ygt(i,:),3,[]);ones(1,nPts)];
    tmp = tmp(1:2,:)./repmat(tmp(3,:),2,1);
    uv(i,:) = reshape(tmp,1,[]);
end

% Add noise to the image locations
rng(12345);
uv = uv + noise*randn(size(uv));

% Compute the groups used for dual decomposition
[G,Eg,W] = compute_groups(3,3,E);

%% Create useful folders and copy files

basehom4ps = 'HOM4PS2/';

for i=1:size(G,1)
    tmp = sprintf('%d',i);
    mkdir(tmp);
    copyfile([basehom4ps 'hom4ps2'],tmp);
    tmp = sprintf('%d/bin',i);
    copyfile([basehom4ps 'bin'],tmp);
end

%% Reconstruct each of the NDefs meshes

fname = 'edge_length_error_group'; % Function that computes the solution for one slave
bname = 'edge_length_error_global'; % Energy function for the whole mesh
plotname = 'compute_and_plot_mesh_los'; % Plot function name

% for def=1:NDefs;

for def=1:100;

    Q = int\[reshape(uv(def,:),2,[]);ones(1,nPts)]; % Lines of sight

    % Ground-truth mesh in camera coordinates
    ygt = ext*[reshape(Ygt(def,:),3,[]);ones(1,nPts)];
        
    % Run DD-Poly or ADMM-Poly
    f = fopen(['./polydata/','model',num2str(def),'.txt'], 'wt+');
    
    fprintf(f, '%6.4f\n', Q);
    fprintf(f, '%6.4f\n', ygt);
    fclose(f);
end
