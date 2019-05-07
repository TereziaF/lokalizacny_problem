function[a, b, c] = presnost(X, pozicie)
q = abs(X - pozicie);
q(1) = q(1)*85.39;                   %prevod stup�ov zem.d�ky na kilometre
q(2) = q(2)*111.03;                  %prevod stup�ov zem.��rky na kilometre
q_of = sqrt(q(1, :).^2 + q(2, :).^2);  %odch�lka
b = max(q_of);                       %maxim�lna odch�lka
a = sum(q_of);                       %s��et absol�tnych odch�lok
c = mean(q_of);                      %priemern� absol�tna odch�lka
end