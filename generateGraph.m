% clc;clear;
% cnputs
% Fully connected graph.
clc;clear;
load vase128.depth
scale = 32;
vase = reshape(vase128, 128, 128);
l = [0, 0, 1];
vase = imresize(vase, 1/scale)/scale;
% vase([1, end],:) = 0;
[vase_x, vase_y] = meshgrid(1:128/scale, 1:128/scale);
surf(vase_x, vase_y, vase); axis equal;
[im, im2, p_gt, q_gt] = imRender(vase_x, vase_y, vase, [0, 0, 1]);
im
im2
imshow(im2);


% im = imresize(im, 1/32);
% vase = imresize(vase/32, 1/32);

[imHeight, imWidth] = size(im);
lambda = 1e-6;
syms x y z c;
fvex = c^4*(5/4*(x-y)^4+5/4*(x-z)^4+2*(x-y)^2*(x-z)^2+2*(x-y)^2+2*(x-z)^2)+(c^2-1)^2 + lambda*(x^2+y^2+z^2);
fcave = -1/4*c^4*((x-y)^4+(x-z)^4)-2*c^2*((x-y)^2+(x-z)^2);
grid_size = imHeight;
vartable = [x y z c];
degree = 2;

hess = 0;

% Generate graph
graph.type = 'triplet';
graph.grid_size = grid_size;
graph.degree = degree;
graph.nodes_num = graph.grid_size^2;
graph.nodes = sym('node', [1 graph.nodes_num]);
if (strcmp(graph.type,'grid'))
    graph.edges_num = 2*(graph.grid_size - 1)*(graph.grid_size);
elseif (strcmp(graph.type,'full'))
    graph.edges_num = graph.nodes_num^2/2 - graph.nodes_num;
elseif (strcmp(graph.type,'triplet'))
    graph.edges_num = (graph.grid_size -2)^2;
else
    graph.edges_num = 1;
end
graph.fpot_sum = sym(0);
graph.fvex_sum = sym(0);
graph.fcave_sum = sym(0);
graph.gpot_sum = sym(0);
graph.gvex_sum = sym(0);
graph.gcave_sum = sym(0);
graph.edges_ind = cell(1, graph.edges_num);
graph.fpot = cell(1, graph.edges_num);
graph.gpot = cell(1, graph.edges_num);
graph.fvex = cell(1, graph.edges_num);
graph.gvex = cell(1, graph.edges_num);
graph.fcave = cell(1, graph.edges_num);
graph.gcave = cell(1, graph.edges_num);
graph.Hpot = cell(1, graph.edges_num);
graph.Hvex = cell(1, graph.edges_num);
graph.Hcave = cell(1, graph.edges_num);

temp_counter = 1;
sub_i = zeros(1,2);
sub_j = zeros(1,2);

for i = 1 : graph.grid_size - 2
    for j = 1 : graph.grid_size - 2
            n1 = sub2ind([grid_size grid_size], i, j);
            n2 = sub2ind([grid_size grid_size], i+1, j);
            n3 = sub2ind([grid_size grid_size], i, j+1);
            fprintf('generate %d th edge...\n', temp_counter);
            graph.edges_ind{temp_counter} = [n1, n2, n3];
            graph.edges{temp_counter} = [graph.nodes(n1),graph.nodes(n2),graph.nodes(n3)];
            graph.fvex{temp_counter} = subs(fvex, vartable, [graph.edges{temp_counter},im(i,j)]);
            graph.fcave{temp_counter} = subs(fcave, vartable, [graph.edges{temp_counter}, im(i,j)]);
            graph.fpot{temp_counter} = graph.fvex{temp_counter}+graph.fcave{temp_counter};
            [graph.gvex{temp_counter}, graph.Hvex{temp_counter}] = cal_hess(graph.fvex{temp_counter}, graph.edges{temp_counter}, hess);
            [graph.gcave{temp_counter}, graph.Hcave{temp_counter}] = cal_hess(graph.fcave{temp_counter}, graph.edges{temp_counter}, hess);
            [graph.gpot{temp_counter}, graph.Hpot{temp_counter}] = cal_hess(graph.fpot{temp_counter}, graph.edges{temp_counter}, hess);
            graph.fpot_sum = graph.fpot_sum + graph.fpot{temp_counter};
            graph.fvex_sum = graph.fvex_sum + graph.fvex{temp_counter};
            graph.fcave_sum = graph.fcave_sum + graph.fcave{temp_counter};
            temp_counter = temp_counter + 1;
    end
end

fprintf('caculate global function...time consuming\n');
[graph.gpot_sum, graph.Hpot_sum] = cal_hess(graph.fpot_sum, graph.nodes, hess);
fprintf('caculate global convex function...time consuming\n');
[graph.gvex_sum, graph.Hvex_sum] = cal_hess(graph.fvex_sum, graph.nodes, hess);
fprintf('caculate global concave function...time consuming\n');
[graph.gcave_sum, graph.Hcave_sum] = cal_hess(graph.fcave_sum, graph.nodes, hess);
fprintf('it is done\n');

save graph16 graph;
