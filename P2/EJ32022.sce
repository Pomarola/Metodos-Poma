clc
clear

function b = miHorner(p, x0)
    n = degree(p) + 1;
    a = coeff(p);
    b = a(n);
    for i = 1:(n-2)
        b = a(n-i) + x0*b;
    end
    b = a(1) + x0*b;
endfunction

function [b,d] = hornerAndDer(p, x0)
    n = degree(p) + 1;
    a = coeff(p);
    b = a(n);
    d = b;
    for i = 1:(n-2)
        b = a(n-i) + x0*b;
        d = b + x0 * d
    end
    b = a(1) + x0*b;
endfunction

p = poly([1 3 2 3 1],"x","coeff");
dp = derivat(p)
[r,q] = hornerAndDer(p, 2);

printf("mi_horner de p en 2: %e\n", r);
printf("horner de p en 2: %e\n", horner(p,2));
printf("mi_horner derivada p: %e\n", q);
printf("horner de dp en 2: %e\n", horner(dp, 2));
