function [func_val, dx] = get_edgeval(z, diff_conv, grad_cave, f_conv, nodes)
func_val = subs(f_conv + grad_cave'*nodes', nodes, z');
dx = subs(diff_conv + grad_cave, nodes, z');
% if (~isLegal(dx))
%     dx = real(dx);
%     % fprintf('error!');
% end
func_val = double(func_val);
dx = double(dx);