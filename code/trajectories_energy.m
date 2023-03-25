%% RANDOM GRAPH
epsilon = 1; %threshold for synchronization
p_edge = 0.2;
N = 20;
adj_matrix = mpm_zeros(N);
for v=1:N
    for w=1:N
        p = rand(1);
        if p < p_edge && v<w
            adj_matrix(v,w) = 0;
            adj_matrix(w,v) = adj_matrix(v,w);
        end
    end
end
G = graph(adj_matrix~=Inf);
M = size(G.Edges,1);
figure(1)
plot(G)
%% RANDOM SHEAF
D = 10; %dimension of matrix weights
a_max = 10;
A = zeros(D,D,N,N);%initialize A
for i=1:N
    for j=1:N
        if adj_matrix(i,j) ~= Inf
            A(:,:,i,j) = a_max*rand(D,D)-a_max;
        else
            A(:,:,i,j) = mp_zeros(D,D);
        end
    end
end
%% HEAT EQUATION
x_max = 10;
p = 1;
%run heat equation
T = 20; % number of iterations
n_trials = 20; %number of trials
E_tarski = zeros(T,1);
iterations = linspace(0,T-1,T);
trace_tarski = zeros(D,N,T,n_trials);
figure; hold on
for trial=1:n_trials
    %initialize X
    X0 = zeros(D,N);
    for i=1:N
        X0(:,i) = x_max*randi(D,1);
    end
    X_tarski = X0;
    E_tarski(1) = dirichlet(A,X_tarski,p);
    trace_tarski(:,:,1,trial)=X_tarski;
    for t=2:T
        X_tarski = mp_add(mpm_add(X_tarski,tarski_laplacian(A,X_tarski,adj_matrix)),mp_ones(D,N));
        trace_tarski(:,:,t,trial) = X_tarski;
        E_tarski(t) = dirichlet(A,X_tarski,p);
    end
    plot(iterations,E_tarski,'black') %plot energy curve
end
plot(iterations,epsilon*M*ones(1,T),'red')



