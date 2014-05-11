% Author: Mathieu Salzmann, NICTA.
%
% Compute the global energy function. This also corresponds to an upper bound on
% the best value we can obtain.
function upb = edge_length_error_global(x,groups,vars)

[Q,E,Eg,W] = deal(vars{:});

f = zeros(1,size(E,1));
for i=1:size(E,1)
    tmp = x(E(i,1))*Q(:,E(i,1)) - x(E(i,2))*Q(:,E(i,2));
    f(i) = tmp'*tmp - E(i,3)^2;
end
upb = 0.5*(f*(f'))/size(E,1);
