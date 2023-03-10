function [Y] = minimax_laplacian(A,X)
    % Tarski Laplacain for assignment X
    % A is the MxDxNxN data array such that A(:,:,i,j) is the mp matrix
    % weight of (directed) edge (i,j)
    % X is a DxN data matrix with D(:,i) a vertex assignment
    [~,D,N,~]=size(A);
    Y = zeros(D,N);
    parfor i=1:N
        Z = mpm_zeros(D,1);
        for j=1:N
            Z = mpm_add(Z,mp_multi(mpm_multi(mp_conv(mp_inv(A(:,:,i,j))),A(:,:,j,i)),X(:,j)));
        end
        Y(:,i) = Z;
    end
end