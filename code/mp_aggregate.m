function [Y] = mean_aggregate(X)
    [Dv,N] = size(X);
    Y = mp_zeros(Dv,1);
    for i=1:N
        Y = mp_add(Y,X(:,i));
    end
end