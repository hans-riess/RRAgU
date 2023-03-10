function [Y] = mpm_aggregate(X)
    [Dv,N] = size(X);
    Y = mpm_zeros(Dv,1);
    for i=1:N
        Y = mpm_add(Y,X(:,i));
    end
end