clc
clear


// Toma un polinomio y un punto
// Devuelta una tupla donde el primer valor es el polinomio evaluado en el punto
// y el segundo la derivada del polinomio evaluado en el mismo punto
function [b,d] = mi_horner(p,x0)
    n = degree(p)+1;
    a = coeff(p)
    b = a(n);
    d = a(n);
    for i=1:(n-2)
        b = a(n-i) + x0*b;
        d = b + x0*d;
    end
    b = a(1) + x0*b;
endfunction

q = poly([1 3 2 3 1], "x", "coeff")
dq = derivat(q)
[r,p] = mi_horner(q, 2);
printf("mi_horner de q en 1: %e\n", r);
printf("horner de q en 1: %e\n", horner(q,2));
printf("mi_horner derivada p: %e\n", p);
printf("horner de q en 1: %e\n", horner(dq, 2));
