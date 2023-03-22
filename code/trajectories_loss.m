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
%% HEAT EQUATION
x_min = 0;
x_max = 10;
D = 3;
p = 1;
%run heat equation
T = 50; % number of iterations
n_trials = 20; %number of trials
E_tarski = zeros(T,1);
iterations = linspace(0,T-1,T);
trace_tarski = zeros(D,N,T,n_trials);
figure; hold on
for trial=1:n_trials
    %initialize X
    X0 = zeros(D,N);
    for i=1:N
        X0(:,i) = randi([x_min,x_max],D,1);
    end
    X_tarski = X0;
    E_tarski(1) = energy_loss(A,X_tarski,p);
    trace_tarski(:,:,1,trial)=X_tarski;
    for t=2:T
        X_tarski = mp_add(mpm_add(X_tarski,tarski_laplacian(A,X_tarski,adj_matrix)),mp_ones(D,N));
        trace_tarski(:,:,t,trial) = X_tarski;
        E_tarski(t) = energy_loss(A,X_tarski,p);
    end
    plot(iterations,E_tarski,'black') %plot loss curve
end
