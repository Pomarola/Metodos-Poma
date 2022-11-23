clc
clear

eps = 0.1

for k=1:10
    mprintf("Matriz %f \n", eps*k)
    A = [1 -1 0; -2 4 -2; 0 -1 (1+eps*k)]
    p = poly(A,"x")
    mprintf("Polinomio caracteristico:")
    disp(p)
    mprintf("Raices del polinomio caracteristico:")
    disp(roots(p))
    printf("Autovalores:")
    disp(spec(A))
end
