function[X, G, Z] = lokalizacia_SDPrelaxacia(N, m, D, Senzory)
n = size(Senzory, 1);         
e = eye(N);

cvx_begin SDP
variable G(N, N) symmetric
variable Z(N+n, N+n) symmetric
variable X(n, N)
minimize 0

%ohranièenie vzdialenosti medzi neznámym a známym/neznámym senzorom
for i = 1:N-m
    for j = 1:N
        if ((i < j) && (D(i,j)~=0))
        H = (e(:,i)-e(:,j))*(e(:,i)-e(:,j))';
        sum(diag(G*H)) == D(i,j);
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

%podmienka semidefinitnosti
Z == semidefinite(n+N);                   
cvx_end

end