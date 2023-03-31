function [Y,aggregates] = tarski_laplacian(A,X,adjacency_matrix)
    % Tarski Laplacain for assignment X
    % A is the MxDxNxN data array such that A(:,:,i,j) is the mp matrix
    % weight of (directed) edge (i,j)
    % X is a DxN data matrix with D(:,i) a vertex assignment
    [~,D,N,~]=size(A);
    Y = zeros(D,N);
    aggregates = zeros(D,N,N);
    parfor i=1:N
        Z = mpm_zeros(D,1);
        for j=1:N
            Z = mpm_add(Z,mpm_multi(adjacency_matrix(i,j),mpm_multi(mp_conv(mp_inv(A(:,:,i,j))), mp_multi(A(:,:,j,i),X(:,j)))));
            aggregates(:,j,i) = mpm_multi(mp_conv(mp_inv(A(:,:,i,j))), mp_multi(A(:,:,j,i),X(:,j)));
        end
        Y(:,i) = Z;
    end
end