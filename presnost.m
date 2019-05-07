function[a, b, c] = presnost(X, pozicie)
q = abs(X - pozicie);
q(1) = q(1)*85.39;                   %prevod stupòov zem.dåžky na kilometre
q(2) = q(2)*111.03;                  %prevod stupòov zem.šírky na kilometre
q_of = sqrt(q(1, :).^2 + q(2, :).^2);  %odchýlka
b = max(q_of);                       %maximálna odchýlka
a = sum(q_of);                       %súèet absolútnych odchýlok
c = mean(q_of);                      %priemerná absolútna odchýlka
end