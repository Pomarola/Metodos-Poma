function [x,a] = gaussElimPPTriDiag(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A (Tridiagonal) y el vector b.
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
    
    if abs(a(i+1,i)) > abs(a(i,i)) then
        temp = a(i,:)
        a(i, :) = a(i+1, :)
        a(i+1, :) = temp
    end
    
    m = a(i+1,i)/a(i,i)
    a(i+1,i) = 0
    a(i+1,i+1) = a(i+1,i+1) - m * a(i,i+1)
    a(i+1,n+1) = a(i+1,n+1) - m * a(i,n+1)
    
end

x(n) = a(n,n+1)/a(n,n);
for i = n-1:-1:1
    x(i) = ( a(i,n+1) - ( a(i,i+1) * x(i+1) ) ) / a(i,i)
end

endfunction

