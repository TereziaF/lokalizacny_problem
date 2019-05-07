%LOKALIZ¡CIA MOBILN›CH TELEF”NOV

%polohy vysielaËov
base_stations = [49.11071 49.11800 49.11156 49.11781 49.11232 49.11919;
20.06544 20.06731 20.06890 20.06859 20.06172 20.06721];

%n·hodnÈ polohy mobilov v okolÌ vysielaËov
mobiles = [(49.12-49.11)*rand(1, 50)+49.11; (20.07-20.06)*rand(1,50)+20.06];

%polohy vöetk˝ch objektov v sieti
pozicie = [mobiles, base_stations];

n = 2;                        %dimenzia priestoru
N = size(pozicie, 2);         %poËet vysielaËov a mobilov
m = size(base_stations, 2);   %poËet vysielaËov

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ZANEDBANIE RUCHU, NEOBMEDZEN› DOSAH SIGN¡LU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%EDM pre vzdialenosti medzi vysielaËmi a mobilmi
D = zeros(N, N);
for i = 1:N
   for j = 1:N     
        D(i, j) = sum((pozicie(:, i)-pozicie(:, j)).^2);
    end
end

%%%%%%%%%%%%%%%%
%SDP relax·cia
%%%%%%%%%%%%%%%%
tic
[X_rel, G_rel, Z_rel] = lokalizacia_SDPrelaxacia(N, m, D, base_stations);
toc
kreslenie(N, m, X_rel, pozicie);
title('SDP relax·cia');
r = hodnost(Z_rel, 10^(-8));
[sucet_odch_rel, max_odch_rel, priem_odch_rel] = presnost(X_rel, pozicie);

disp(['SDP relax·cia: hodnosù ', num2str(r)]) 
disp(['SDP relax·cia: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_rel)])
disp(['SDP relax·cia: maxim·lna odch˝lka(km) ', num2str(max_odch_rel)])
disp(['SDP relax·cia: priemern· odch˝lka(km) ', num2str(priem_odch_rel)])

%%%%%%%%%%%%%%%%%%%
%Heuristika stopy
%%%%%%%%%%%%%%%%%%%
tic
[X_stopa, G_stopa, Z_stopa] = lokalizacia_heuristika_stopy(N, m, D, base_stations);
toc
kreslenie(N, m, X_stopa, pozicie);
title('Heuristika stopy');
r = hodnost(Z_stopa, 10^(-8));
[sucet_odch_stopa, max_odch_stopa, priem_odch_stopa] = presnost(X_stopa, pozicie);

disp(['Heuristika stopy: hodnosù ', num2str(r)]) 
disp(['Heuristika stopy: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_stopa)])
disp(['Heuristika stopy: maxim·lna odch˝lka(km) ', num2str(max_odch_stopa)])
disp(['Heuristika stopy: priemern· odch˝lka(km) ', num2str(priem_odch_stopa)])

%%%%%%%%%%%%%%%%%%%%
%Konvexn· iter·cia
%%%%%%%%%%%%%%%%%%%%
tic
%prv· iter·cia - heuristika stopy
[X_kvx, G_kvx, Z_kvx] = lokalizacia_heuristika_stopy(N, m, D, base_stations);
%prv· iter·cia - SDP relax·cia
%[X_kvx, G_kvx, Z_kvx] = lokalizacia_SDP_relaxacia(N, m, D, base_stations);
W = lokalizacia_konvexna_iteracia2(Z_kvx, N, n);
iter=1;

while ((hodnost(Z_kvx, 10^(-8))~=2) && (iter<20))
[X_kvx, G_kvx, Z_kvx] = lokalizacia_konvexna_iteracia1(N, m, D, base_stations, W);
W = lokalizacia_konvexna_iteracia2(Z_kvx, N, n);
iter=iter+1;
end
toc

r = hodnost(Z_kvx, 10^(-8));
kreslenie(N, m,  X_kvx, pozicie);
title('Konvexn· iter·cia');
[sucet_odch_kvx, max_odch_kvx, priem_odch_kvx] = presnost(X_kvx, pozicie);


disp(['Konvexn· iter·cia: hodnosù ', num2str(r)]) 
disp(['Konvexn· iter·cia: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_kvx)])
disp(['Konvexn· iter·cia: maxim·lna odch˝lka(km) ', num2str(max_odch_kvx)])
disp(['Konvexn· iter·cia: priemern· odch˝lka(km) ', num2str(priem_odch_kvx)])
disp(['Konvexn· iter·cia: PoËet iter·ciÌ(km) ', num2str(iter)]) 

%%%%%%%%%%%%%%%%
%Heuristika 1a
%%%%%%%%%%%%%%%%

hodnost_1a = zeros(1, 9);
sucet_odch_1a = zeros(1, 9);
max_odch_1a = zeros(1, 9);
priem_odch_1a = zeros(1, 9);

%voæba optim·lnej relatÌvnej v·hy
for beta = 0:0.25:2
[X_1a,G_1a,Z_1a] = lokalizacia_heuristika_1a(N, m, D, base_stations, beta);
hodnost_1a(int16(beta*4+1)) = hodnost(Z_1a, 10^(-8));
[sucet_odch_1a(int16(beta*4+1)), max_odch_1a(int16(beta*4+1)), priem_odch_1a(int16(beta*4+1))] = presnost(X_1a, pozicie);
end

koeficienty = [0:0.25:2];
beta = koeficienty(find(max_odch_1a == min(max_odch_1a)));

tic
[X_1a,G_1a,Z_1a] = lokalizacia_heuristika_1a(N, m, D, base_stations, beta);
toc
r = hodnost(Z_1a, 10^(-8));
[sucet_odch_1a, max_odch_1a, priem_odch_1a] = presnost(X_1a, pozicie);
kreslenie(N, m, X_1a, pozicie);
title('Heuristika 1a');

disp(['Heuristika 1a: hodnosù ', num2str(r)]) 
disp(['Heuristika 1a: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_1a)])
disp(['Heuristika 1a: maxim·lna odch˝lka(km) ', num2str(max_odch_1a)])
disp(['Heuristika 1a: priemern· odch˝lka(km) ', num2str(priem_odch_1a)])

%%%%%%%%%%%%%%%%
%Heuristika 1b
%%%%%%%%%%%%%%%%
tic
[X_1b, G_1b, Z_1b] = lokalizacia_heuristika_1b(N, m, D, base_stations);
toc
r = hodnost(Z_1b, 10^(-8));
[sucet_odch_1b, max_odch_1b, priem_odch_1b] = presnost(X_1b, pozicie);
kreslenie(N, m, X_1b, pozicie);
title('Heuristika 1b');

disp(['Heuristika 1b: hodnosù ', num2str(r)]) 
disp(['Heuristika 1b: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_1b)])
disp(['Heuristika 1b: maxim·lna odch˝lka(km) ', num2str(max_odch_1b)])
disp(['Heuristika 1b: priemern· odch˝lka(km) ', num2str(priem_odch_1b)])

%%%%%%%%%%%%%%%%
%Heuristika 2a
%%%%%%%%%%%%%%%%

hodnost_2a = zeros(1, 9);
sucet_odch_2a = zeros(1, 9);
max_odch_2a = zeros(1, 9);
priem_odch_2a = zeros(1, 9);

%voæba optim·lnej relatÌvnej v·hy
for beta = 0:0.25:2
[X_2a,G_2a,Z_2a] = lokalizacia_heuristika_2a(N, m, D, base_stations, beta);
hodnost_2a(int16(beta*4+1)) = hodnost(Z_2a, 10^(-8));
[sucet_odch_2a(int16(beta*4+1)), max_odch_2a(int16(beta*4+1)), priem_odch_2a(int16(beta*4+1))] = presnost(X_2a, pozicie);
end

koeficienty = [0:0.25:2];
beta = koeficienty(find(max_odch_2a == min(max_odch_2a)));
tic
[X_2a, G_2a, Z_2a] = lokalizacia_heuristika_2a(N, m, D, base_stations, beta);
toc
r = hodnost(Z_2a, 10^(-8));
[sucet_odch_2a, max_odch_2a, priem_odch_2a] = presnost(X_2a, pozicie);
kreslenie(N, m, X_2a, pozicie);
title('Heuristika 2a');

disp(['Heuristika 2a: hodnosù ', num2str(r)]) 
disp(['Heuristika 2a: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_2a)])
disp(['Heuristika 2a: maxim·lna odch˝lka(km) ', num2str(max_odch_2a)])
disp(['Heuristika 2a: priemern· odch˝lka(km) ', num2str(priem_odch_2a)])

%%%%%%%%%%%%%%%%
%Heuristika 2b
%%%%%%%%%%%%%%%%
tic
[X_2b, G_2b, Z_2b] = lokalizacia_heuristika_2b(N, m, D, base_stations);
toc
r = hodnost(Z_2b, 10^(-8));
[sucet_odch_2b, max_odch_2b, priem_odch_2b] = presnost(X_2b, pozicie);
kreslenie(N, m, X_2b, pozicie);
title('Heuristika 2b');

disp(['Heuristika 2b: hodnosù ', num2str(r)]) 
disp(['Heuristika 2b: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_2b)])
disp(['Heuristika 2b: maxim·lna odch˝lka(km) ', num2str(max_odch_2b)])
disp(['Heuristika 2b: priemern· odch˝lka(km) ', num2str(priem_odch_2b)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ZANEDBANIE RUCHU, OBMEDZEN› DOSAH SIGN¡LU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%voæba dosahu sign·lu (v stupÚoch zem.öÌrky)
R = 0.005;
%R = 0.0075;
 
%EDM pre vzdialenosti medzi vysielaËmi a mobilmi
A = zeros(N, N);
for i = 1:N
    for j = 1:N 
        A(i, j) = sum((pozicie(:, i)-pozicie(:, j)).^2);
    end
end
D = zeros(N, N);
for i = 1:N
    for j = 1:N 
        if(A(i, j) <= R^2 || (i > N-m && j > N-m))
        D(i, j) = sum((pozicie(:, i)-pozicie(:, j)).^2);
        end
    end
end

%percento zn·mych prvkov EDM matice
percento = (length(D(D>0))/2)/(N*(N-1)/2);

%%%%%%%%%%%%%%%%
%SDP relax·cia
%%%%%%%%%%%%%%%%
tic
[X_rel, G_rel, Z_rel] = lokalizacia_SDPrelaxacia(N, m, D, base_stations);
toc
kreslenie(N, m, X_rel, pozicie);
title('SDP relax·cia');
r = hodnost(Z_rel, 10^(-8));
[sucet_odch_rel, max_odch_rel, priem_odch_rel] = presnost(X_rel, pozicie);

disp(['SDP relax·cia: hodnosù ', num2str(r)]) 
disp(['SDP relax·cia: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_rel)])
disp(['SDP relax·cia: maxim·lna odch˝lka(km) ', num2str(max_odch_rel)])
disp(['SDP relax·cia: priemern· odch˝lka(km) ', num2str(priem_odch_rel)])

%%%%%%%%%%%%%%%%%%%
%Heuristika stopy
%%%%%%%%%%%%%%%%%%%
tic
[X_stopa, G_stopa, Z_stopa] = lokalizacia_heuristika_stopy(N, m, D, base_stations);
toc
kreslenie(N, m, X_stopa, pozicie);
title('Heuristika stopy');
r = hodnost(Z_stopa, 10^(-8));
[sucet_odch_stopa, max_odch_stopa, priem_odch_stopa] = presnost(X_stopa, pozicie);

disp(['Heuristika stopy: hodnosù ', num2str(r)]) 
disp(['Heuristika stopy: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_stopa)])
disp(['Heuristika stopy: maxim·lna odch˝lka(km) ', num2str(max_odch_stopa)])
disp(['Heuristika stopy: priemern· odch˝lka(km) ', num2str(priem_odch_stopa)])

%%%%%%%%%%%%%%%%%%%%
%Konvexn· iter·cia
%%%%%%%%%%%%%%%%%%%%
tic
%prv· iter·cia - heuristika stopy
[X_kvx, G_kvx, Z_kvx] = lokalizacia_heuristika_stopy(N, m, D, base_stations);
%prv· iter·cia - SDP relax·cia
%[X_kvx, G_kvx, Z_kvx] = lokalizacia_SDP_relaxacia(N, m, D, base_stations);
W = lokalizacia_konvexna_iteracia2(Z_kvx, N, n);
iter=1;

while ((hodnost(Z_kvx, 10^(-8))~=2) && (iter<20))
[X_kvx, G_kvx, Z_kvx] = lokalizacia_konvexna_iteracia1(N, m, D, base_stations, W);
W = lokalizacia_konvexna_iteracia2(Z_kvx, N, n);
iter=iter+1;
end
toc

r = hodnost(Z_kvx, 10^(-8));
kreslenie(N, m, X_kvx, pozicie);
title('Konvexn· iter·cia');
[sucet_odch_kvx, max_odch_kvx, priem_odch_kvx] = presnost(X_kvx, pozicie);


disp(['Konvexn· iter·cia: hodnosù ', num2str(r)]) 
disp(['Konvexn· iter·cia: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_kvx)])
disp(['Konvexn· iter·cia: maxim·lna odch˝lka(km) ', num2str(max_odch_kvx)])
disp(['Konvexn· iter·cia: priemern· odch˝lka(km) ', num2str(priem_odch_kvx)])
disp(['Konvexn· iter·cia: PoËet iter·ciÌ(km) ', num2str(iter)]) 

%%%%%%%%%%%%%%%%
%Heuristika 1a
%%%%%%%%%%%%%%%%

hodnost_1a = zeros(1, 9);
sucet_odch_1a = zeros(1, 9);
max_odch_1a = zeros(1, 9);
priem_odch_1a = zeros(1, 9);

%voæba optim·lnej relatÌvnej v·hy
for beta = 0:0.25:2
[X_1a,G_1a,Z_1a] = lokalizacia_heuristika_1a(N, m, D, base_stations, beta);
hodnost_1a(int16(beta*4+1)) = hodnost(Z_1a, 10^(-8));
[sucet_odch_1a(int16(beta*4+1)), max_odch_1a(int16(beta*4+1)), priem_odch_1a(int16(beta*4+1))] = presnost(X_1a, pozicie);
end

koeficienty = [0:0.25:2];
beta = koeficienty(find(max_odch_1a == min(max_odch_1a)));

tic
[X_1a, G_1a, Z_1a] = lokalizacia_heuristika_1a(N, m, D, base_stations, beta);
toc
r = hodnost(Z_1a, 10^(-8));
[sucet_odch_1a, max_odch_1a, priem_odch_1a] = presnost(X_1a, pozicie);
kreslenie(N, m, X_1a, pozicie);
title('Heuristika 1a');

disp(['Heuristika 1a: hodnosù ', num2str(r)]) 
disp(['Heuristika 1a: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_1a)])
disp(['Heuristika 1a: maxim·lna odch˝lka(km) ', num2str(max_odch_1a)])
disp(['Heuristika 1a: priemern· odch˝lka(km) ', num2str(priem_odch_1a)])

%%%%%%%%%%%%%%%%
%Heuristika 1b
%%%%%%%%%%%%%%%%
tic
[X_1b, G_1b, Z_1b] = lokalizacia_heuristika_1b(N, m, D, base_stations);
toc
r = hodnost(Z_1b, 10^(-8));
[sucet_odch_1b, max_odch_1b, priem_odch_1b] = presnost(X_1b, pozicie);
kreslenie(N, m, X_1b, pozicie);
title('Heuristika 1b');

disp(['Heuristika 1b: hodnosù ', num2str(r)]) 
disp(['Heuristika 1b: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_1b)])
disp(['Heuristika 1b: maxim·lna odch˝lka(km) ', num2str(max_odch_1b)])
disp(['Heuristika 1b: priemern· odch˝lka(km) ', num2str(priem_odch_1b)])

%%%%%%%%%%%%%%%%
%Heuristika 2a
%%%%%%%%%%%%%%%%

hodnost_2a = zeros(1, 9);
sucet_odch_2a = zeros(1, 9);
max_odch_2a = zeros(1, 9);
priem_odch_2a = zeros(1, 9);

%voæba optim·lnej relatÌvnej v·hy
for beta = 0:0.25:2
[X_2a,G_2a,Z_2a] = lokalizacia_heuristika_2a(N, m, D, base_stations, beta);
hodnost_2a(int16(beta*4+1)) = hodnost(Z_2a, 10^(-8));
[sucet_odch_2a(int16(beta*4+1)), max_odch_2a(int16(beta*4+1)), priem_odch_2a(int16(beta*4+1))] = presnost(X_2a, pozicie);
end

koeficienty = [0:0.25:2];
beta = koeficienty(find(max_odch_2a == min(max_odch_2a)));
tic
[X_2a,G_2a,Z_2a] = lokalizacia_heuristika_2a(N, m, D, base_stations, beta);
toc
r = hodnost(Z_2a, 10^(-8));
[sucet_odch_2a, max_odch_2a, priem_odch_2a] = presnost(X_2a, pozicie);
kreslenie(N, m, X_2a, pozicie);
title('Heuristika 2a');

disp(['Heuristika 2a: hodnosù ', num2str(r)]) 
disp(['Heuristika 2a: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_2a)])
disp(['Heuristika 2a: maxim·lna odch˝lka(km) ', num2str(max_odch_2a)])
disp(['Heuristika 2a: priemern· odch˝lka(km) ', num2str(priem_odch_2a)])

%%%%%%%%%%%%%%%%
%Heuristika 2b
%%%%%%%%%%%%%%%%
tic
[X_2b, G_2b, Z_2b] = lokalizacia_heuristika_2b(N, m, D, base_stations);
toc
r = hodnost(Z_2b, 10^(-8));
[sucet_odch_2b, max_odch_2b, priem_odch_2b] = presnost(X_2b, pozicie);
kreslenie(N, m, X_2b, pozicie);
title('Heuristika 2b');

disp(['Heuristika 2b: hodnosù ', num2str(r)]) 
disp(['Heuristika 2b: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_2b)])
disp(['Heuristika 2b: maxim·lna odch˝lka(km) ', num2str(max_odch_2b)])
disp(['Heuristika 2b: priemern· odch˝lka(km) ', num2str(priem_odch_2b)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRÕTOMNOSç RUCHU, NEOBEMEDZEN› DOSAH SIGN¡LU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ni=0.001;              %ruchov˝ faktor
D_lower = zeros(N, N);  %EDM matica doln˝ch odhadov vzdialenostÌ
D_upper = zeros(N, N);  %EDM matica horn˝ch odhadov vzdialenostÌ

%simul·cia doln˝ch a horn˝ch odhadov vzdialenostÌ
for i = 1:N
    %vzdialenosti medzi mobilmi a mobilmi/vysielaËmi s˙ ovplyvnenÈ ruchom
    for j = 1:N
    D_lower(i, j) = sum((pozicie(:, i)-pozicie(:, j)).^2)*(1-ni*rand(1));
    D_upper(i, j) = sum((pozicie(:, i)-pozicie(:, j)).^2)*(1+ni*rand(1));
    end
end
    %vzdialenosti medzi vysielaËmi s˙ presnÈ
for i = N-m+1:N
    for j = N-m+1:N 
        D_lower(i, j) = sum((pozicie(:, i)-pozicie(:, j)).^2); 
        D_upper(i, j) = D_lower(i, j);
    end
end

%%%%%%%%%%%%%%%%
%SDP relax·cia
%%%%%%%%%%%%%%%%
tic
[X_rel, G_rel, Z_rel] = lokalizacia_SDPrelaxacia_ruch(N, m, D_lower, D_upper, base_stations);
toc
kreslenie(N, m, X_rel, pozicie);
title('SDP relax·cia, ruch');
r = hodnost(Z_rel, 10^(-8));
[sucet_odch_rel, max_odch_rel, priem_odch_rel] = presnost(X_rel,pozicie);

disp(['SDP relax·cia: hodnosù ', num2str(r)]) 
disp(['SDP relax·cia: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_rel)])
disp(['SDP relax·cia: maxim·lna odch˝lka(km) ', num2str(max_odch_rel)])
disp(['SDP relax·cia: priemern· odch˝lka(km) ', num2str(priem_odch_rel)])

%%%%%%%%%%%%%%%%%%%
%Heuristika stopy
%%%%%%%%%%%%%%%%%%%
tic
[X_stopa, G_stopa, Z_stopa] = lokalizacia_heuristika_stopy_ruch(N, m, D_lower, D_upper, base_stations);
toc
r = hodnost(Z_stopa, 10^(-8));
kreslenie(N, m, X_stopa, pozicie);
title('Heuristika stopy, ruch');
[sucet_odch_stopa, max_odch_stopa, pocet_odch_stopa] = presnost(X_stopa, pozicie);

disp(['Heuristika stopy: hodnosù ', num2str(r)]) 
disp(['Heuristika stopy: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_stopa)])
disp(['Heuristika stopy: maxim·lna odch˝lka(km) ', num2str(max_odch_stopa)])
disp(['Heuristika stopy: priemern· odch˝lka(km) ', num2str(priem_odch_stopa)])

%%%%%%%%%%%%%%%%%%%%
%Konvexn· iter·cia
%%%%%%%%%%%%%%%%%%%%
tic
[X_kvx, G_kvx, Z_kvx] = lokalizacia_heuristika_stopy_ruch(N, m, D_lower, D_upper, base_stations);
W = lokalizacia_konvexna_iteracia2(Z_kvx, N, n);
iter=1;

while ((hodnost(Z_kvx, 10^(-8))~=2) && (iter<20))
[X_kvx,G_kvx, Z_kvx] = lokalizacia_konvexna_iteracia1_ruch(N, m, D_lower, D_upper, base_stations, W);
U = lokalizacia_konvexna_iteracia2(Z_kvx, N, n);
iter=iter+1;
end
toc

r = hodnost(Z_kvx, 10^(-8));
kreslenie(N, m, X_kvx, pozicie);
title('Konvexn· iter·cia, ruch')
[sucet_odch_kvx, max_odch_kvx, priem_odch_kvx] = presnost(X_kvx,pozicie);

disp(['Konvexn· iter·cia: hodnosù ', num2str(r)]) 
disp(['Konvexn· iter·cia: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_kvx)])
disp(['Konvexn· iter·cia: maxim·lna odch˝lka(km) ', num2str(max_odch_kvx)])
disp(['Konvexn· iter·cia: priemern· odch˝lka(km) ', num2str(priem_odch_kvx)])
disp(['Konvexn· iter·cia: PoËet iter·ciÌ(km) ', num2str(iter)]) 

%%%%%%%%%%%%%%%%
%Heuristika 1a
%%%%%%%%%%%%%%%%

hodnost_1a = zeros(1, 9);
sucet_odch_1a = zeros(1, 9);
max_odch_1a = zeros(1, 9);
priem_odch_1a = zeros(1, 9);

%voæba optim·lnej relatÌvnej v·hy
for beta = 0:0.25:2
[X_1a,G_1a,Z_1a] = lokalizacia_heuristika_1a_ruch(N, m, D_lower, D_upper, base_stations, beta);
hodnost_1a(int16(beta*4+1)) = hodnost(Z_1a, 10^(-8));
[sucet_odch_1a(int16(beta*4+1)), max_odch_1a(int16(beta*4+1)), priem_odch_1a(int16(beta*4+1))] = presnost(X_1a, pozicie);
end

koeficienty = [0:0.25:2];
beta = koeficienty(find(max_odch_1a == min(max_odch_1a)));
tic
[X_1a,G_1a,Z_1a] = lokalizacia_heuristika_1a_ruch(N, m, D_lower, D_upper, base_stations, beta);
toc

r = hodnost(Z_1a, 10^(-8));
[sucet_odch_1a, max_odch_1a, priem_odch_1a] = presnost(X_1a, pozicie);
kreslenie(N, m, X_1a, pozicie);
title('Heuristika 1a, ruch');

disp(['Heuristika 1a: hodnosù ', num2str(r)]) 
disp(['Heuristika 1a: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_1a)])
disp(['Heuristika 1a: maxim·lna odch˝lka(km) ', num2str(max_odch_1a)])
disp(['Heuristika 1a: priemern· odch˝lka(km) ', num2str(priem_odch_1a)])

%%%%%%%%%%%%%%%%
%Heuristika 1b
%%%%%%%%%%%%%%%%
tic
[X_1b, G_1b, Z_1b] = lokalizacia_heuristika_1b_ruch(N, m, D_lower, D_upper, base_stations);
toc

r = hodnost(Z_1b, 10^(-8));
[sucet_odch_1b, max_odch_1b, pocet_odch_1b] = presnost(X_1b, pozicie);
kreslenie(N, m, X_1b, pozicie);
title('Heuristika 1b, ruch');

disp(['Heuristika 1b: hodnosù ', num2str(r)]) 
disp(['Heuristika 1b: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_1b)])
disp(['Heuristika 1b: maxim·lna odch˝lka(km) ', num2str(max_odch_1b)])
disp(['Heuristika 1b: priemern· odch˝lka(km) ', num2str(priem_odch_1b)])

%%%%%%%%%%%%%%%%
%heuristika 2a
%%%%%%%%%%%%%%%%

hodnost_2a = zeros(1, 9);
sucet_odch_2a = zeros(1, 9);
max_odch_2a = zeros(1, 9);
priem_odch_2a = zeros(1, 9);

%voæba optim·lnej relatÌvnej v·hy
for beta = 0:0.25:2
[X_2a,G_2a,Z_2a] = lokalizacia_heuristika_2a_ruch(N, m, D_lower, D_upper, base_stations, beta);
hodnost_2a(int16(beta*4+1)) = hodnost(Z_2a, 10^(-8));
[sucet_odch_2a(int16(beta*4+1)), max_odch_2a(int16(beta*4+1)), priem_odch_2a(int16(beta*4+1))] = presnost(X_2a, pozicie);
end

koeficienty = [0:0.25:2];
beta = koeficienty(find(max_odch_2a == min(max_odch_2a)));
tic
[X_2a, G_2a, Z_2a] = lokalizacia_heuristika_2a_ruch(N, m, D_lower, D_upper, base_stations, beta);
toc
r = hodnost(Z_2a, 10^(-8));
[sucet_odch_2a, max_odch_2a, priem_odch_2a] = presnost(X_2a, pozicie);
kreslenie(N, m, X_2a, pozicie);
title('Heuristika 2a, ruch')

disp(['Heuristika 2a: hodnosù ', num2str(r)]) 
disp(['Heuristika 2a: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_2a)])
disp(['Heuristika 2a: maxim·lna odch˝lka(km) ', num2str(max_odch_2a)])
disp(['Heuristika 2a: priemern· odch˝lka(km) ', num2str(priem_odch_2a)])

%%%%%%%%%%%%%%%%
%Heuristika 2b
%%%%%%%%%%%%%%%%
tic
[X_2b, G_2b, Z_2b] = lokalizacia_heuristika_2b_ruch(N, m, D_lower, D_upper, base_stations);
toc
r = hodnost(Z_2b, 10^(-8));
kreslenie(N, m, X_2b, pozicie);
title('Heuristika 2b, ruch');
[sucet_odch_2b, max_odch_2b, priem_odch_2b] = presnost(X_2b, pozicie);

disp(['Heuristika 2b: hodnosù ', num2str(r)]) 
disp(['Heuristika 2b: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_2b)])
disp(['Heuristika 2b: maxim·lna odch˝lka(km) ', num2str(max_odch_2b)])
disp(['Heuristika 2b: priemern· odch˝lka(km) ', num2str(priem_odch_2b)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRÕTOMNOSç RUCHU, OBMEDZEN› DOSAH SIGN¡LU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%voæba dosahu sign·lu (v stupÚoch zem.öÌrky)
R = 0.005;
%R=0.0075;
 
%EDM pre vzdialenosti medzi vysielaËmi a mobilmi
C = zeros(N, N);
for i = 1:N
    for j = 1:N 
        C(i, j) = sum((pozicie(:, i)-pozicie(:, j)).^2);
    end
end
D_lower_new = zeros(N, N);
D_upper_new = zeros(N, N);
for i = 1:N
    for j = 1:N 
        if(C(i, j) <= R^2 || (i > N-m && j > N-m))
        D_lower_new(i, j) = D_lower(i, j);
        D_upper_new(i, j) = D_upper(i, j);
        end
    end
end

%%%%%%%%%%%%%%%%
%SDP relax·cia
%%%%%%%%%%%%%%%%
tic
[X_rel, G_rel, Z_rel] = lokalizacia_SDPrelaxacia_ruch(N, m, D_lower_new, D_upper_new, base_stations);
toc
kreslenie(N, m, X_rel, pozicie);
title('SDP relax·cia, ruch');
r = hodnost(Z_rel, 10^(-8));
[sucet_odch_rel, max_odch_rel, priem_odch_rel] = presnost(X_rel, pozicie);

disp(['SDP relax·cia: hodnosù ', num2str(r)]) 
disp(['SDP relax·cia: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_rel)])
disp(['SDP relax·cia: maxim·lna odch˝lka(km) ', num2str(max_odch_rel)])
disp(['SDP relax·cia: priemern· odch˝lka(km) ', num2str(priem_odch_rel)])

%%%%%%%%%%%%%%%%%%%
%Heuristika stopy
%%%%%%%%%%%%%%%%%%%
tic
[X_stopa, G_stopa, Z_stopa] = lokalizacia_heuristika_stopy_ruch(N, m, D_lower_new, D_upper_new, base_stations);
toc
r = hodnost(Z_stopa, 10^(-8));
kreslenie(N, m, X_stopa, pozicie);
title('Heuristika stopy, ruch');
[sucet_odch_stopa, max_odch_stopa, priem_odch_stopa] = presnost(X_stopa, pozicie);

disp(['Heuristika stopy: hodnosù ', num2str(r)]) 
disp(['Heuristika stopy: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_stopa)])
disp(['Heuristika stopy: maxim·lna odch˝lka(km) ', num2str(max_odch_stopa)])
disp(['Heuristika stopy: priemern· odch˝lka(km) ', num2str(priem_odch_stopa)])

%%%%%%%%%%%%%%%%%%%%
 

%%%%%%%%%%%%%%%%
%Heuristika 1a
%%%%%%%%%%%%%%%%

hodnost_1a = zeros(1, 9);
sucet_odch_1a = zeros(1, 9);
max_odch_1a = zeros(1, 9);
priem_odch_1a = zeros(1, 9);

%voæba optim·lnej relatÌvnej v·hy
for beta = 0:0.25:2
[X_1a, G_1a, Z_1a] = lokalizacia_heuristika_1a_ruch(N, m, D_lower_new, D_upper_new, base_stations, beta);
hodnost_1a(int16(beta*4+1)) = hodnost(Z_1a, 10^(-8));
[sucet_odch_1a(int16(beta*4+1)), max_odch_1a(int16(beta*4+1)), priem_odch_1a(int16(beta*4+1))] = presnost(X_1a, pozicie);
end

koeficienty = [0:0.25:2];
beta = koeficienty(find(max_odch_1a == min(max_odch_1a)));
tic
[X_1a, G_1a, Z_1a] = lokalizacia_heuristika_1a_ruch(N, m, D_lower_new, D_upper_new, base_stations, beta);
toc

r = hodnost(Z_1a, 10^(-8));
[sucet_odch_1a, max_odch_1a, priem_odch_1a] = presnost(X_1a, pozicie);
kreslenie(N, m, X_1a, pozicie);
title('Heuristika 1a, ruch');

disp(['Heuristika 1a: hodnosù ', num2str(r)]) 
disp(['Heuristika 1a: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_1a)])
disp(['Heuristika 1a: maxim·lna odch˝lka(km) ', num2str(max_odch_1a)])
disp(['Heuristika 1a: priemern· odch˝lka(km) ', num2str(priem_odch_1a)])

%%%%%%%%%%%%%%%%
%Heuristika 1b
%%%%%%%%%%%%%%%%
tic
[X_1b, G_1b, Z_1b] = lokalizacia_heuristika_1b_ruch(N, m, D_lower_new, D_upper_new, base_stations);
toc

r = hodnost(Z_1b, 10^(-8));
[sucet_odch_1b, max_odch_1b, pocet_odch_1b] = presnost(X_1b, pozicie);
kreslenie(N, m, X_1b, pozicie);
title('Heuristika 1b, ruch');

disp(['Heuristika 1b: hodnosù ', num2str(r)]) 
disp(['Heuristika 1b: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_1b)])
disp(['Heuristika 1b: maxim·lna odch˝lka(km) ', num2str(max_odch_1b)])
disp(['Heuristika 1b: priemern· odch˝lka(km) ', num2str(priem_odch_1b)])

%%%%%%%%%%%%%%%%
%heuristika 2a
%%%%%%%%%%%%%%%%

hodnost_2a = zeros(1, 9);
sucet_odch_2a = zeros(1, 9);
max_odch_2a = zeros(1, 9);
priem_odch_2a = zeros(1, 9);

%voæba optim·lnej relatÌvnej v·hy
for beta = 0:0.25:2
[X_2a, G_2a, Z_2a] = lokalizacia_heuristika_2a_ruch(N, m, D_lower_new, D_upper_new, base_stations, beta);
hodnost_2a(int16(beta*4+1)) = hodnost(Z_2a, 10^(-8));
[sucet_odch_2a(int16(beta*4+1)), max_odch_2a(int16(beta*4+1)), priem_odch_2a(int16(beta*4+1))] = presnost(X_2a, pozicie);
end

koeficienty = [0:0.25:2];
beta = koeficienty(find(max_odch_2a == min(max_odch_2a)));
tic
[X_2a,G_2a,Z_2a] = lokalizacia_heuristika_2a_ruch(N, m, D_lower_new, D_upper_new, base_stations, beta);
toc
r = hodnost(Z_2a, 10^(-8));
[sucet_odch_2a, max_odch_2a, priem_odch_2a] = presnost(X_2a, pozicie);
kreslenie(N, m, X_2a, pozicie);
title('Heuristika 2a, ruch')

disp(['Heuristika 2a: hodnosù ', num2str(r)]) 
disp(['Heuristika 2a: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_2a)])
disp(['Heuristika 2a: maxim·lna odch˝lka(km) ', num2str(max_odch_2a)])
disp(['Heuristika 2a: priemern· odch˝lka(km) ', num2str(priem_odch_2a)])

%%%%%%%%%%%%%%%%
%Heuristika 2b
%%%%%%%%%%%%%%%%
tic
[X_2b, G_2b, Z_2b] = lokalizacia_heuristika_2b_ruch(N, m, D_lower_new, D_upper_new, base_stations);
toc
r = hodnost(Z_2b, 10^(-8));
kreslenie(N, m, X_2b, pozicie);
title('Heuristika 2b, ruch');
[sucet_odch_2b, max_odch_2b, priem_odch_2b] = presnost(X_2b,pozicie);

disp(['Heuristika 2b: hodnosù ', num2str(r)]) 
disp(['Heuristika 2b: s˙Ëet odch˝lok(km) ', num2str(sucet_odch_2b)])
disp(['Heuristika 2b: maxim·lna odch˝lka(km) ', num2str(max_odch_2b)])
disp(['Heuristika 2b: priemern· odch˝lka(km) ', num2str(priem_odch_2b)])