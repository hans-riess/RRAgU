%% RANDOM GRAPH
p_edge = 0.5;
N = 10;
adj_matrix = mpm_zeros(N);
for v=1:N
    for w=1:N
        p = rand(1);
        if p > p_edge && v<w
            adj_matrix(v,w) = floor(5*(1-p));
            adj_matrix(w,v) = floor(5*(p-p));
        end
    end
end
G = graph(adj_matrix~=Inf);
figure(1)
plot(G)
%% RANDOM SHEAF
Dv = 3;
De = 3;
a_min = -10;
a_max = 1;
%initialize A
A = zeros(De,Dv,N,N);
for i=1:N
    for j=1:N
        if adj_matrix(i,j) ~= Inf
            A(:,:,i,j) = randi([a_min,a_max],De,Dv);
        else
            A(:,:,i,j) = mp_zeros(De,Dv);
        end
    end
end
%% INITIAL CONDITION
x_min = 0;
x_max = 10;
Dv = 3;
%initialize X
X0 = zeros(Dv,N);
for i=1:N
    X0(:,i) = randi([x_min,x_max],Dv,1);
end
%% HEAT EQUATION
p = Inf;
%run heat equation
T = 50; % number of iterations
E_tarski = zeros(T,1);
X_tarski = X0;
trace_tarski = zeros(Dv,N,T);
E_tarski(1) = dirichlet(A,X_tarski,p);
trace_tarski(:,:,1)=X_tarski;
for t=2:T
    X_tarski = mpm_add(X_tarski,tarski_laplacian(A,X_tarski,adj_matrix));
    trace_tarski(:,:,t) = X_tarski;
    E_tarski(t) = dirichlet(A,X_tarski,p);
end

%plot the results
figure(2)
plot(linspace(0,T-1,T),E_tarski)
