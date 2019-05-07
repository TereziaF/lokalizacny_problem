function[X, G, Z] = lokalizacia_SDPrelaxacia(N, m, D, Senzory)
n = size(Senzory, 1);         
e = eye(N);

cvx_begin SDP
variable G(N, N) symmetric
variable Z(N+n, N+n) symmetric
variable X(n, N)
minimize 0

%ohrani�enie vzdialenosti medzi nezn�mym a zn�mym/nezn�mym senzorom
for i = 1:N-m
    for j = 1:N
        if ((i < j) && (D(i,j)~=0))
        H = (e(:,i)-e(:,j))*(e(:,i)-e(:,j))';
        sum(diag(G*H)) == D(i,j);
        end
    end
end

%ohrani�enie d�ok zn�mych senzorov
for k = N-m+1:N
  F = e(:,k)*e(:,k)';
  sum(diag(G*F)) == Senzory(:,k-N+m)'*Senzory(:,k-N+m);  
end

%ohrani�enie vzdialenosti medzi zn�mymi senzormi
for i = N-m+1:N
    for j = N-m+1:N
        J = (e(:,i)*e(:,j)' + e(:,j)*e(:,i)')/2;
        if i < j
            sum(diag(G*J)) == Senzory(:,i-N+m)'*Senzory(:,j-N+m);
        end
    end
end

%zadanie zn�mych senzorov
X(:, N-m+1:N) == Senzory;   

%defin�cia matice Z
Z == [eye(n), X; X', G];

%podmienka semidefinitnosti
Z == semidefinite(n+N);                   
cvx_end

end