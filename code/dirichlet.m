function [energy] = dirichlet(A,X,p)
    [De,Dv,N,~] = size(A);
    energy = 0;
    for i=1:N
        for j=1:N
            if any(any(A(:,:,i,j)> mp_zeros(De,Dv)))
                energy = energy + norm(mp_multi(A(:,:,j,i),X(:,j))-mp_multi(A(:,:,i,j),X(:,i)),p);
            end
        end
    end
    energy = 1/2*energy;
end