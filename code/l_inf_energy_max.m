function [energy] = l_inf_energy_max(adjacency_matrix,p)
    [N,~] = size(adjacency_matrix);
    energy = 0;
    for i=1:N
        for j=1:N
            if adjacency_matrix(i,j)~=Inf
                energy = energy + adjacency_matrix(i,j);
            end
        end
    end
    energy = 0.5*energy;
end
