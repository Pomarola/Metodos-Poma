clear
clc

function [x,a] = gaussElimMultiple(A,B)
// Esta función obtiene la solución del sistema de multiples ecuaciones lineales A*x=B, 
// dada la matriz de coeficientes A y la matriz de soluciones.
// La función implementa el método de Eliminación Gaussiana sin pivoteo.  

[nA,mA] = size(A) 
[nB,mB] = size(B)

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nB then
    error('gausselim - dimensiones incompatibles entre A y B');
    abort;
end;

a = [A B]; // Matriz aumentada
bStart = nA+1
bEnd = nA+mB

// Eliminación progresiva

for i = 1:nA-1
    for j= i+1:nA
        mji = a(j,i)/a(i,i)
        a(j,i) = 0
        a(j,i+1:bEnd) = a(j,i+1:bEnd) - mji * a(i,i+1:bEnd)
    end
end

x(nA,:) = a(nA,bStart:bEnd)/a(nA,nA);
for i = nA-1:-1:1
    x(i,:) = ( a(i,bStart:bEnd) - ( a(i, i+1:nA) * x(i+1:nA,:) ) ) / a(i,i)
end

endfunction


// Ejemplos de aplicación

A = [1 2 3; 3 -2 1; 4 2 -1]
b1 = [14 2 5]'
b2 = [9 -5 19]'
b3 = [-2 2 12]'
B = [b1 b2 b3]

// x =
//   1.   2.   2.
//   2.   5.   1.
//   3.  -1.  -2.
[x,a] = gaussElimMultiple(A,B)
disp(x)
disp(a)

A = [1 2 3; 3 -2 1; 4 2 -1]
b1 = [1 0 0]'
b2 = [0 1 0]'
b3 = [0 0 1]'
B = [b1 b2 b3]

// x =
//   0.      0.1428571   0.1428571
//   0.125  -0.2321429   0.1428571
//   0.25    0.1071429  -0.1428571
[x,a] = gaussElimMultiple(A,B)
disp(x)
disp(a)
