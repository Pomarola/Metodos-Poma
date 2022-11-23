function [x,a] = gaussElimPP(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana con pivoteo.  

[nA,mA] = size(A) 
[nb,mb] = size(b)

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;

a = [A b]; // Matriz aumentada

// Eliminación progresiva
n = nA;
for i = 1:n-1
    pivot = i
    amax = abs(a(i,i))  //pivoteo
    for j=i+1:n
        if abs(a(j,i)) > amax then
            pivot = j
            amax = a(j,i)
        end
    end
    temp = a(pivot,:)
    a(pivot,:) = a(i,:)
    a(i,:) = temp
    
    for j=i+1:n
        mji = a(j,i)/a(i,i)
        a(j,i) = 0
        a(j,i+1:n+1) = a(j,i+1:n+1) - mji * a(i,i+1:n+1)
    end
end

x(n) = a(n,n+1)/a(n,n);
for i = n-1:-1:1
    x(i) = ( a(i,n+1) - ( a(i, i+1:n) * x(i+1:n) ) ) / a(i,i)
end

endfunction


// Ejemplos de aplicación

A = [1 1 0 3; 2 1 -1 1; 3 -1 -1 2; -1 2 3 -1]
b = [4 1 -3 4]'

// x = [-1 2 0.00000 1]'
[x,a] = gaussElimPP(A,b)
disp(x)
disp(a)

A2 = [1 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
b2 = [-8 -20 -2 4]'

// x = [-7 3 2 2]'
[x2,a2] = gaussElimPP(A2,b2)
disp(x2)
disp(a2)

A3 = [1 1 0 4; 2 1 -1 1; 4 -1 -2 2; 3 -1 -1 2]
b3 = [2 1 0 -3]'

// x = [-4 2/3 -7 4/3]'
[x3,a3] = gaussElimPP(A3,b3)
disp(x3)
disp(a3)
