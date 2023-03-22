%% RANDOM NUMBER GENERATOR
% rng(1)

%% RANDOM GRAPH
D = 3; %dimension of matrix weights
epsilon = 10; %threshold for synchronization
p_edge = 0.5;
N = 10;
adj_matrix = mpm_zeros(N);
for v=1:N
    for w=1:N
        p = rand(1);
        if p < p_edge && v<w
            adj_matrix(v,w) = epsilon/D+(2*rand(1)-1);
            adj_matrix(w,v) = epsilon/D+(2*rand(1)-1);
        end
    end
end
G = graph(adj_matrix~=Inf);
figure(1)
plot(G)
%% RANDOM SHEAF
a_min = -10;
a_max = 1;
%initialize A
A = zeros(D,D,N,N);
for i=1:N
    for j=1:N
        if adj_matrix(i,j) ~= Inf
            A(:,:,i,j) = randi([a_min,a_max],D,D);
        else
            A(:,:,i,j) = mp_zeros(D,D);
        end
    end
end
%% CENTRALIZED SOLUTION


