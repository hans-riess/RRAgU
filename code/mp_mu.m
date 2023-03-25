function [mu] = mp_mu(A,b)
    mu = 0.5*norm(mp_multi(A,mpm_multi(mp_conv(mp_inv(A)),b))-b,Inf);
end