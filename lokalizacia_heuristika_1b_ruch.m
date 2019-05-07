function[X, G, Z] = lokalizacia_heuristika_1b_ruch(N, m, D_lower, D_upper, Senzory)
n = size(Senzory, 1);       
e = eye(N);
pocet = n;

cvx_begin SDP
variable G(N, N) symmetric
variable Z(N+n, N+n) symmetric
variable Y(2*N+2*n, 2*N+2*n) symmetric
variable X(n, N)
variable z(1, 1)
variable V(N+n, N+n) symmetric
minimize pocet*z + sum(diag(V))

%ohranièenie vzdialenosti medzi neznámym a známym/neznámym senzorom
for i = 1:N-m
    for j = 1:N
        if ((i < j)&&(D_lower(i,j)~=0)&&(D_upper(i,j)~=0))
        H = (e(:,i)-e(:,j))*(e(:,i)-e(:,j))';
        D_lower(i,j) <= sum(diag(G*H)) <= D_upper(i,j);
        end
    end
end

%ohranièenie dåžok známych senzorov
for k = N-m+1:N
  F = e(:,k)*e(:,k)';
  sum(diag(G*F)) == Senzory(:,k-N+m)'*Senzory(:,k-N+m);  
end

%ohranièenie vzdialenosti medzi známymi senzormi
for i = N-m+1:N
    for j = N-m+1:N
        J = (e(:,i)*e(:,j)' + e(:,j)*e(:,i)')/2;
        if i < j
            sum(diag(G*J)) == Senzory(:,i-N+m)'*Senzory(:,j-N+m);
        end
    end
end

%zadanie známych senzorov
X(:, N-m+1:N) == Senzory; 

%definícia matice Z
Z == [eye(n), X; X', G];

%podmienka z SDP pre súèet vl.hodnôt
z*eye(n+N, n+N) + V - Z == semidefinite(n+N);

%nové ohranièenie - súèet n najväèších vl.hodnôt rovný stope
pocet*z + sum(diag(V)) == sum(diag(Z));

%definícia matice Y
Y == [Z zeros(N+n, N+n); zeros(N+n, N+n) V];

%podmienka semidefinitnosti
Y == semidefinite(2*n+2*N);  

cvx_end
end