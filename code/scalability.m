%% EXPERIMENT
%start pool
poolStartup

%random number generator
rng(1)

%parameters
x_min = -10;
x_max = 10;
N_min = 10;
N_max = 1000;
N_trials = 20;
NN = ceil(linspace(N_min,N_max,N_trials));
x_min = -10;
x_max = 10;
epsilon = 5; %threshold for synchronization
p_edge = 0.2;
D_min = 2;
D_max = 20;
D_trials = 5;
DD = ceil(linspace(D_min,D_max,D_trials));
scalable_matrix = zeros(N_trials,D_trials);
%% Experiment
for trial_N=19:N_trials
    for trial_D=4:D_trials
        N = ceil(NN(trial_N));
        D = ceil(DD(trial_D));
        %generate graph
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
        %generate sheaf
        A = zeros(D,D,N,N);
        for i=1:N
            for j=1:N
                if adj_matrix(i,j) ~= Inf
                    A(:,:,i,j) = randi([x_min,x_max],D,D);
                else
                    A(:,:,i,j) = mp_zeros(D,D);
                end
            end
        end
        %generate initial state
        X = zeros(D,N);
        for i=1:N
            X(:,i) = randi([x_min,x_max],D,1);
        end
        %run Tarski Laplacin & time
        f = @() tarski_laplacian(A,X,adj_matrix);
        scalable_matrix(trial_N,trial_D) = timeit(f);
    end
end
%% PLOT RESULTS
labels = 'D = ' + string(reshape(DD,D_trials,1));
figure,hold on;
title('Scalability Analysis')
xlabel('Number of Nodes')
ylabel('Time (s)')
for trial_D=1:D_trials
    plot(NN(1:18),scalable_matrix(1:18,trial_D),'--s')
end
legend(labels)
