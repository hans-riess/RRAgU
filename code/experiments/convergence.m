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
M = size(G.Edges,1);
%% RANDOM SHEAF
D = 10; %dimension of matrix weights
A = zeros(D,D,N,N);%initialize A
for agent=1:N
    for j=1:N
        if adj_matrix(agent,j) ~= Inf
            A(:,:,agent,j) = 2*rand(D,D)-1;
        else
            A(:,:,agent,j) = mp_zeros(D,D);
        end
    end
end
%% HEAT EQUATION
p = Inf;
%run heat equation
T = 10; % number of iterations
n_trials = 20; %number of trials
loss_tarski = zeros(T+1,n_trials);
iterations = linspace(0,T,T+1);
trace_tarski = zeros(D,N,T+1,n_trials);
for trial=1:n_trials
    %initialize X
    X0 = zeros(D,N);
    for agent=1:N
        X0(:,agent) = 2*rand(D,1)-1;
    end
    X_tarski = X0;
    loss_tarski(1,trial) = energy_loss(A,X_tarski,p);
    trace_tarski(:,:,1,trial)=X_tarski;
    for t=2:T+1
        X_tarski = mp_add(mpm_add(X_tarski,tarski_laplacian(A,X_tarski,adj_matrix)),mp_ones(D,N));
        trace_tarski(:,:,t,trial) = X_tarski;
        loss_tarski(t,trial) = energy_loss(A,X_tarski,p);
    end
end

%% CONVERGENCE ANALYSIS

alpha = zeros(T,n_trials);
X = reshape(trace_tarski,[N*D,T+1,n_trials]);
for trial=1:n_trials
    for t=1:T
        alpha(t,trial) = norm(X(:,t+1,trial)-X(:,t,trial),"inf");
    end
end

%% PLOTS

close all

%colors
orange = "#ff7f00";
blue = "#377eb8";
red = "#e41a1c";
green = "#4daf4a";

%PLOT1
non_zero_alpha = [3,7,17,20];
figure; hold on
for trial=1:n_trials
    if ismember(trial,non_zero_alpha)
        plot(iterations(1:T),alpha(:,trial),"Color",orange,'LineWidth',1)
    else
        plot(iterations(1:T),alpha(:,trial),"Color",blue,'LineWidth',1)
    end
end
title('Convergence')
xlabel('Iterations ($t$)','Interpreter','latex')
ylabel('$\alpha(t)$','Interpreter','latex')
axis([0 9 -0.05 1])
yticks([0])
legend('$\alpha=0$','','$\alpha>0$','Interpreter','latex')

%PLOT 2
figure; hold on
plot(iterations,loss_max(adj_matrix)*ones(T+1),'Color',red,'LineStyle','--','LineWidth',1.5)
for trial=1:n_trials
    if ismember(trial,non_zero_alpha)
        plot(iterations,loss_tarski(:,trial),'Color',orange,"LineWidth",1)
    else
        plot(iterations,loss_tarski(:,trial),'Color',blue,"LineWidth",1)
    end
end
title('Loss')
xlabel('Iterations ($t$)','Interpreter','latex')
yticks([loss_max(adj_matrix)])
ylabel('$\ell\bigl(\mathbf{X}(t)\bigr)$','Interpreter','latex')
legend('$\epsilon$','Interpreter','latex')

%% AGENTS

close all
trial = 3; %trial with alpha<0
trace = trace_tarski(:,:,:,trial); %trace of that trial
for agent=1:N
    figure; hold on
    title('Agent ' + string(agent))
    for alternative=1:D
        values = reshape(trace(alternative,agent,:),[T+1,1]);
        plot(iterations,values,'Color',orange)
    end
end

%% GRAPH

figure
bad_nodes = [4,9,14];
h = plot(G,"Layout",'force','NodeColor','black','EdgeColor','black');
path = shortestpath(G,4,14);
highlight(h,bad_nodes)
highlight(h,path,'NodeColor',orange,'EdgeColor',orange,'LineWidth',1.5)
title('Non-converging Agents')


