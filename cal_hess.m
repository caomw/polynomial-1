function [f_grad, f_hess] = cal_hess(f1, vartable, opt)

f_grad = vartable';
for i = 1 : length(vartable)
    f_grad(i) = [diff(f1, vartable(i))];
end

if (opt)
    f_hess = repmat(f_grad, [1 length(vartable)]);
    for i = 1 : length(vartable)
        for j = 1 : length(vartable)
            f_hess(i,j) = diff(f_grad(i),vartable(j));
        end
    end
else
    f_hess = 0;
end
