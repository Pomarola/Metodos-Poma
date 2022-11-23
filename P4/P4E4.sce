clear
clc

function d = gaussDet(A)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana sin pivoteo.  

[nA,mA] = size(A) 

if nA<>mA then
    d = 0
    return;
end

// Eliminación progresiva
n = nA;
for i = 1:n-1
    for j= i+1:n
        mji = A(j,i)/A(i,i)
        A(j,i) = 0
        A(j,i+1:n) = A(j,i+1:n) - mji * A(i,i+1:n)
    end
end

d = 1;
for i = 1:n
    d = d*A(i,i)
end

endfunction


// Ejemplos de aplicación

A = [1 1 0 3; 2 1 -1 1; 3 -1 -1 2; -1 2 3 -1]

// d = 39
d = gaussDet(A)
disp(d)

A2 = [5 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]

// d = -116
d = gaussDet(A2)
disp(d)

A3 = [1 0 0; 0 1 0; 0 0 1]
// d = 1
d = gaussDet(A3)
disp(d)


A4 = [1 0 0 1; 0 1 1 0; 0 0 1 1]
// d = 0
d = gaussDet(A4)
disp(d)
