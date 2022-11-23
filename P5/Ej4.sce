clc
clear 

function [A,b] = create4(N)
    A = 8*eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1) + diag(ones(N-3,1),3) + diag(ones(N-3,1),-3)
    b = ones(N,1)
endfunction

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

//Metodo de GaussSeidel iterativo
//Recibe la matriz A, el vector solucion b, un vector inicial x0, y la tolerancia eps
function x = gaussSeidel(A, b, x0, eps)
    n = size(A, 1)
    x = x0
    xk = x0
    cont = 0
    suma = 0
    
    for i=1:n-1
        [v,index]=max(abs(A(i:n,i)))
        pivot = i-1+index
        temp = A(pivot,:); A(pivot,:) = A(i,:)
        A(i,:) = temp
        temp = b(pivot,:); b(pivot,:) = b(i,:)
        b(i,:) = temp
    end
    
    I = eye(n,n)
    N = tril(A)
    N = I-inv(N)*A
    rspec = max(abs(spec(N)))
    if (rspec < 1) then
        disp("Converge")
    else
        disp("No Converge")
        x = %nan
        return
    end
    
    while (norm(x-xk) > eps || cont == 0)
        xk = x
        for i = 1:n
            suma = 0
            for j = 1:i-1
                suma = suma + A(i,j)*x(j)
            end
            for j = i+1:n
                suma = suma + A(i,j)*x(j)
            end
            x(i) = 1/A(i,i) * (b(i)-suma)
        end
        cont = cont+1
    end
    disp(cont)
endfunction

[A,b] = create4(100)
tic()
[a, x] = gaussElimPP(A,b)
t=toc()
disp(t) // 0.111159
//disp(x)

tic()
x = gaussSeidel(A,b,zeros(100,1),0.000001)
t=toc()
disp(t)
// Converge
//
//   18.
//
//   0.7545484
//disp(x)

tic()
x = gaussSeidel(A,b,zeros(100,1),0.00000000001)
t=toc()
disp(t)
// Converge
//
//   40.
//
//   1.5046611
//disp(x)

[A,b] = create4(500)
tic()
[a, x] = gaussElimPP(A,b)
t=toc()
disp(t) // 4.1443442
//disp(x)

tic()
x = gaussSeidel(A,b,zeros(500,1),0.000001)
t=toc()
disp(t)
// Converge
//
//   18.
//
//   16.545329
//disp(x)

tic()
x = gaussSeidel(A,b,zeros(500,1),0.00000000001)
t=toc()
disp(t)
// Converge
//
//   40.
//
//   36.701728
//disp(x)
