function a = kreslenie(N, m, X, pozicie)

%vykreslenie skuto�n�ch poz�ci�
for i = 1:N
    a = circle(pozicie(1, i), pozicie(2, i), 5*10^(-4));
    hold on
end
xlabel('Zemepisn� d�ka')
ylabel('Zemepisn� ��rka')

%vykreslenie n�jden�ch poz�ci�
for i = N-m+1:N
       hold on
       a = plot(X(1, i), X(2, i), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
end
for i = 1:N-m
       hold on
       a = plot(X(1, i), X(2, i), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
end
axis equal
set(gca, 'box', 'off')

end