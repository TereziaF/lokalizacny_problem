function[W] = lokalizacia_konvexna_iteracia2(X, N, n)

cvx_begin SDP
variable W(n+N, n+N)
minimize sum(diag(W*X))
sum(diag(W)) == N;
W == semidefinite(N+n);
eye(N+n, N+n)-W == semidefinite(N+n);
cvx_end

end