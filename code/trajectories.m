%% RANDOM GRAPH
p_edge = 0.2;
N = 20;
adj_matrix = mpm_zeros(N);
for v=1:N
    for w=1:N
        p = rand(1);
        if p < p_edge && v<w
            adj_matrix(v,w) = rand(1);
            adj_matrix(w,v) = adj_matrix(v,w);
        end
    end
end
G = graph(adj_matrix~=Inf);
M = size(G.Edges,1);
%% RANDOM SHEAF
D = 10; %dimension of matrix weights
A = zeros(D,D,N,N);%initialize A
for i=1:N
    for j=1:N
        if adj_matrix(i,j) ~= Inf
            A(:,:,i,j) = 2*rand(D,D)-1;
        else
            A(:,:,i,j) = mp_zeros(D,D);
        end
    end
end
%% HEAT EQUATION
p = Inf;
%run heat equation
T = 10; % number of iterations
n_trials = 20; %number of trials
E_tarski = zeros(T+1,n_trials);
loss_tarski = zeros(T+1,n_trials);
E_tarski_adapt = zeros(T+1,n_trials);
iterations = linspace(0,T,T+1);
trace_tarski = zeros(D,N,T+1,n_trials);
trace_tarski_adapt = zeros(D,N,T+1,n_trials);
for trial=1:n_trials
    %initialize X
    X0 = zeros(D,N);
    for i=1:N
        X0(:,i) = 2*rand(D,1)-1;
    end
    X_tarski = X0;
    X_tarski_adapt = X0;
    E_tarski_adapt(1,trial) = dirichlet(A,X_tarski_adapt,p);
    E_tarski(1,trial) = dirichlet(A,X_tarski,p);
    loss_tarski(1,trial) = energy_loss(A,X_tarski,p);
    trace_tarski_adapt(:,:,1,trial)=X_tarski;
    trace_tarski(:,:,1,trial)=X_tarski;
    for t=2:T+1
        X_tarski = mp_add(mpm_add(X_tarski,tarski_laplacian(A,X_tarski,adj_matrix)),mp_ones(D,N));
        X_tarski_adapt = mp_add(mpm_add(X_tarski_adapt,tarski_laplacian_adapt(A,X_tarski_adapt)),mp_ones(D,N));
        trace_tarski(:,:,t,trial) = X_tarski;
        E_tarski(t,trial) = dirichlet(A,X_tarski,p);
        loss_tarski(t,trial) = energy_loss(A,X_tarski,p);
        trace_tarski_adapt(:,:,t,trial) = X_tarski_adapt;
        E_tarski_adapt(t,trial) = dirichlet(A,X_tarski_adapt,p);
    end
end
%% PLOT 1
figure; hold on
for trial=1:n_trials
%     plot(iterations,E_tarski_adapt(:,trial),'green') %plot energy curve
    plot(iterations,E_tarski(:,trial),'green')
%     plot(iterations,l_inf_energy_max(adj_matrix)*ones(size(iterations)),'red--');
end
title('Dirichlet energy of trajectories of the heat equation')
xlabel('Iterations')
ylabel('Dirichlet energy')
% legend('Adaptive')
% legend('Weighted')
%% PLOT 2
figure; hold on
plot(iterations,get_epsilon(adj_matrix)*ones(T+1),'red--')
for trial=1:n_trials
    plot(iterations,loss_tarski(:,trial),'black')
end
title('Loss of trajectories of the heat equation')
xlabel('Iterations')
ylabel('Loss')
legend('$\epsilon$','Interpreter','latex')

%% CONVERGENCE ANALYSIS

alpha = zeros(T,n_trials);
X = reshape(trace_tarski,[N*D,T+1,n_trials]);
for trial=1:n_trials
    for t=1:T
        alpha(t,trial) = norm(X(:,t+1,trial)-X(:,t,trial),"inf");
    end
end

%% PLOT 3
non_zero = [3,7,17,20];
figure; hold on
for trial=1:n_trials
    if ismember(trial,non_zero)
        plot(iterations(1:T),alpha(:,trial),"red")
    else
        plot(iterations(1:T),alpha(:,trial),"black")
    end
end
title('Convergence Analysis')
xlabel('$t$','Interpreter','latex')
ylabel('$\| \mathbf{X}(t+1) - \mathbf{X}(t) \|_{\infty}$','Interpreter','latex')
axis([0 9 -0.05 1])
yticks([0 1])
legend('$\alpha=0$','','$\alpha>0$','Interpreter','latex')
%% PLOT 4
figure
plot(G)



