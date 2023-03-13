%% RANDOM SHEAF

%seed
% seed(10)

% G is a ring graph with 3 nodes and 3 edges
G = [0,1,1;1,0,1;1,1,0];
N = 3;
Dv = 3;
De = 3;

%initialize X
X0 = zeros(Dv,N);
for i=1:N
    X0(:,i) = 10*rand(Dv,1);
end

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
T = 20; % number of iterations
E_tarski = zeros(T,1);
X_tarski = X0;
trace_tarski = zeros(Dv,N,T);
E_tarski(1) = dirichlet(A,X_tarski,p);
trace_tarski(:,:,1)=X_tarski;
for t=2:T
    X_tarski = mpm_add(X_tarski,tarski_laplacian(A,X_tarski));
    trace_tarski(:,:,t) = X_tarski;
    E_tarski(t) = dirichlet(A,X_tarski,p);
end

%plot the results
plot(linspace(0,T-1,T),E_tarski)
