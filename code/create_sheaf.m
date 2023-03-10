%% RANDOM SHEAF

%seed
% seed(10)

% G is a ring graph with 3 nodes and 3 edges
G = [0,1,1;1,0,1;1,1,0];
N = 3;
Dv = 3;
De = 3;

%initialize X
a_min = 0;
a_max = 10;
X_tarski = zeros(Dv,N);
for i=1:N
    X_tarski(:,i) = randi([a_min,a_max],Dv,1);
end

X_tarski = [1,7,7;1,5,4;8,3,5];

%initialize A
A = zeros(De,Dv,N,N);
for i=1:N
    for j=1:N
        if G(i,j) == 1
            A(:,:,i,j) = randi([a_min,a_max],De,Dv);
        elseif G(i,j) == 0
            A(:,:,i,j) = mp_zeros(De,Dv);
        end
    end
end
%% HEAT EQUATION
p = Inf;
%run heat equation
T = 5g0; % number of iterations
E_tarski = zeros(T,1);
E_minimax = zeros(T,1);
X_tarski = X0;
X_minimax = X0;
trace_tarski = zeros(Dv,N,T);
trace_minimax = zeros(Dv,N,T);
E_tarski(1) = dirichlet(A,X_tarski,p);
E_minimax(1) = dirichlet(A,X_tarski,p);
trace_tarski(:,:,1)=X_tarski;
trace_minimax(:,:,1) = X_minimax;
for t=2:T
    X_tarski = mpm_add(X_tarski,tarski_laplacian(A,X_tarski));
    X_minimax = mpm_add(X_minimax,minimax_laplacian(A,X_minimax));
    trace_tarski(:,:,t) = X_tarski;
    trace_minimax(:,:,t) = X_minimax;
    E_tarski(t) = dirichlet(A,X_tarski,p);
    E_minimax(t) = dirichlet(A,X_minimax,p);
end

%plot the results
plot(linspace(0,T-1,T),E_tarski,linspace(0,T-1,T),E_minimax)
