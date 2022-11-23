clc
clear

//Metodo de GaussSeidel SOR iterativo
//Recibe la matriz A, el vector solucion b, un vector inicial x0, la tolerancia eps y un factor de relajacion w
function x = SOR(A, b, x0, eps, w)
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
            x(i) = (1-w)*x(i) + w/A(i,i) * (b(i)-suma)
        end
        cont = cont+1
    end
    disp(cont)
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

//Dada una matriz tridiagonal definida positiva, obtiene el parametro optimo para el metodo SOR
function w = optTriDiag(A)
    n = size(A,1)
    I = eye(n,n)
    N = diag(1./diag(A))
    N = I-N*A
    rspec = max(abs(spec(N)))
    w = 2/(1 + sqrt(1 - rspec ^ 2))
endfunction

A = [4 3 0; 3 4 -1; 0 -1 4]
b = [24 30 -24]'
x0 = [0 0 0]'
//x = [3 4 -5]'
eps1 = 10^-2
eps2 = 10^-7
x = SOR(A, b, x0, eps1, optTriDiag(A))
disp(x)
// Converge
//
//   7.
//
//   3.0023348
//   3.9989237
//  -5.0004475

x = gaussSeidel(A, b, x0, eps2)
disp(x)
// Converge
//
//   36.
//
//   3.0000001
//   3.9999999
//  -5.

x = SOR(A, b, x0, eps2, optTriDiag(A))
disp(x)
// Converge
//
//   16.
//
//   3.
//   4.
//  -5.



