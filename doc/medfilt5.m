%%
x = rand(1000,1);
d1 = [x(1); x(1:end-1)];
d2 = [d1(1); d1(1:end-1)];
d3 = [d2(1); d2(1:end-1)];
d4 = [d3(1); d3(1:end-1)];
plot(x); hold on;
plot(d4);
hold off;

A = max(x, d1);  B = min(x, d1);
C = max(d2, d3); D = min(d2, d3);
E = max(B, D);   D = min(B, D);
B = max(C, d4);  C = min(C, d4);
F = max(A, B);   B = min(A, B);
A = max(E, C);   C = min(E, C);
E = max(A, B);   B = min(A, B);
D = max(D, C);  %C = min(D, C);
B = max(B, D);

truth = medfilt1(x, 5);
med_err = truth(1:end-2) - B(3:end);
plot(x); hold on; 
plot(truth(1:end-2)); plot(B(3:end)); hold off;
plot(med_err);