% function [nodes_val, f_sum] = cccp_polynomial_graph_lbfgs(graph, nodes_val, flag_sol)

% flag_sol = 1; % 1 for Mathieu, 0 for cccp
flag_dual = 0;
options.Method = 'lbfgs';
options.Display = 'on';
options.MaxIter = 10;
tic;
nodes_num = graph.nodes_num;
edges_num = graph.edges_num;
% nodes_val = rand(1, nodes_num);
% nodes_val = reshape(vase+rand(size(vase)), 1, nodes_num);
% nodes_val(1) = 0;
k = ceil((size(vase, 1)+1) / 2);
[x, y] = meshgrid(1:size(vase, 1), 1:size(vase, 1));
nodes_val = (k - sqrt((x-(size(vase, 1)+1) / 2).^2+(y-(size(vase, 1)+1)/2).^2))/k*2;

nodes_old = randn(1, nodes_num);
gvex_sum = graph.gvex_sum;
gcave_sum = graph.gcave_sum;
gpot_sum = graph.gpot_sum;

lambda = cell(1, edges_num);
weights = zeros(1, nodes_num);
for i = 1 : edges_num
    lambda{i} = zeros(1, 2);
    weights(graph.edges_ind{i}) = weights(graph.edges_ind{i}) + 1;
end

max_iter = 200;
f_sum = zeros(1, max_iter);
f_old = 0;

for iter = 1 : max_iter
    f_sum(iter) = double(subs(graph.fpot_sum, graph.nodes, nodes_val));
    fprintf('iter %d, current energy %f: \n ', iter, f_sum(iter));
    t = 0;
    nodes_iter = zeros(1, nodes_num);
    if flag_dual
        for i = 1 : edges_num
            z_init = nodes_val(graph.edges_ind{i});
            grad_cave = subs(graph.gcave{i}, graph.edges{i}, z_init);
            z = minFunc(@get_edgeval,z_init',options, graph.gvex{i}, grad_cave, graph.fvex{i}, graph.edges{i});
            z = real(z);
            nodes_iter(graph.edges_ind{i}) = nodes_iter(graph.edges_ind{i}) + z';
            fprintf('.');
        end
        nodes_val = real(nodes_iter ./ (weights+eps));
    else
        z_init = nodes_val;
        grad_cave = subs(graph.gcave_sum, graph.nodes, z_init);
        z = minFunc(@get_edgeval,z_init',options, graph.gvex_sum, grad_cave, graph.fvex_sum, graph.nodes);
        nodes_val = real(z)';
    end
    if (norm(nodes_old - nodes_val) < 1e-8 || abs(f_sum(iter) - f_old) < 1e-8)
        break;
    end
    % reshape(nodes_val, 8, 8)
    f_old = f_sum(iter);
    nodes_old = nodes_val;
    save temp nodes_val;
end
fprintf('\n');
toc
