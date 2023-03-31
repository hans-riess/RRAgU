function l = loss_max(adjacency_matrix)
[N,~] = size(adjacency_matrix);
    l = 0;
    for i=1:N
        for j=1:N
            if adjacency_matrix(i,j)~=Inf
                l = max(l,adjacency_matrix(i,j));
            end
        end
    end
end