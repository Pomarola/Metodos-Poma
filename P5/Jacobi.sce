clc
clear

//Metodo de Jacobi iterativo
//Recibe la matriz A, el vector solucion b, un vector inicial x0, y la tolerancia eps
function x = jacobi(A, b, x0, eps)
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
    N = diag(diag(A))
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
                suma = suma + A(i,j)*xk(j)
            end
            for j = i+1:n
                suma = suma + A(i,j)*xk(j)
            end
            x(i) = 1/A(i,i) * (b(i)-suma)
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

// Ej1
A = [0 2 4; 1 -1 -1; 1 -1 2]
b = [0 0.375 0]'
x0 = [0 0 0]'
x = jacobi(A, b, x0, 0.01)
disp(x)
// Solucion esperada
// No Converge
//
//   Nan

A = [1 -1 0; -1 2 -1; 0 -1 1.1]
b = [0 1 0]'
x0 = [0 0 0]'
x = jacobi(A, b, x0, 0.01)
disp(x)
// Solucion esperada
// Converge
//   171. iteraciones esperadas
//
//   10.789086 solucion esperada
//   10.798673
//   9.8082602

A = [0 2 4; 1 -1 -1; 1 -1 2]
b = [0 0.375 0]'
x0 = [0 0 0]'
x = gaussSeidel(A, b, x0, 0.01)
disp(x)
// Solucion esperada
// No Converge
//
//   Nan

A = [1 -1 0; -1 2 -1; 0 -1 1.1]
b = [0 1 0]'
x0 = [0 0 0]'
x = gaussSeidel(A, b, x0, 0.01)
disp(x)
// Solucion esperada
// Converge
//
//   97. iteraciones esperadas
//
//   10.873565 solucion esperada
//   10.879312
//   9.8902836

// Ej2
A = [10 1 2 3 4; 1 9 -1 2 -3; 2 -1 7 3 -5; 3 2 3 12 -1; 4 -3 -5 -1 15]
b = [12 -27 14 -17 12]'
x0 = [0 0 0 0 0]'
x = jacobi(A, b, x0, 0.000001)
disp(x)
// Solucion esperada
// Converge
// 
//   67. iteraciones esperadas
//
//   1.0000016 solucion esperada
//  -2.0000015
//   2.9999973
//  -1.9999996
//   0.9999981

x = gaussSeidel(A, b, x0, 0.000001)
disp(x)
// Solucion esperada
// Converge
//
//   38. iteraciones esperadas
//
//   1.0000009 solucion esperada
//  -2.0000007
//   2.9999987
//  -1.9999999
//   0.9999992
