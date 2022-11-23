function x = regresive_sust(A, b)
    n = size(A,1)
    a = [A b]
    x(n) = a(n,n+1)/a(n,n);
    for i = n-1:-1:1
        x(i) = ( a(i,n+1) - ( a(i,i+1) * x(i+1) ) ) / a(i,i)
    end     
endfunction

function sol = progresive_sust(A, b)
    n = size(A,1)
    x(1) = b(1)/A(1,1);
    for i = 2:n
        x(i) = ( b(i) - ( A(i, 1:i-1) * x(1:i-1) ) ) / A(i,i)
    end
    sol = x
endfunction



// Metodo Doolittle de factorizacion LU
function [L,U]= doolittle(A)
    [nA, mA] = size(A)
    if nA<>mA then
        error('doolittle - La matriz A debe ser cuadrada');
        abort;
    end
    
    L = eye(A)
    U = eye(A)
    for i=1:nA
        for j=i:nA
            suma = 0
            for k=1:i-1
                suma = suma + L(i,k)*U(k,j)
            end
            U(i,j)=A(i,j) - suma
            for m=i+1:nA
                suma = 0
                for k=1:i-1
                    suma = suma + L(m,k)*U(k,i)
                end
                L(m,i) = (A(m,i)-suma)/U(i,i)
            end
        end
    end
endfunction

// Resuelve un sistema de ecuaciones aplicando Doolittle
function [L,U,x]= d_solver(A, b)
    [L,U] = doolittle(A)
    y = progresive_sust(L,b)
    x = regresive_sust(U,y)
endfunction
// FALTA METODO DE CROUT

function [L,U] = crout(A)
    [nA, mA] = size(A)
    if nA<>mA then
        error('crout - La matriz A debe ser cuadrada');
        abort;
    end
    
    L = eye(A)
    U = eye(A)
    for i=1:nA
        for j=i:nA
            suma = 0
            for k=1:i-1
                suma = suma + L(j,k)*U(k,i)
            end
            L(j,i)=A(j,i) - suma
            for m=i+1:nA
                suma = 0
                for k=1:i-1
                    suma = suma + L(i,k)*U(k,m)
                end
                U(i,m) = (A(i,m)-suma)/L(i,i)
            end
        end
    end
endfunction

function [L,U,x] = crout_solver(A,b)
    [L,U]=crout(A)
    y = progresive_sust(L,b)
    x = regresive_sust(U,y)
endfunction

A = [1 2 3 4;1 4 9 16;1 8 27 64; 1 16 81 256]
b = [2 ;10; 44; 190]

[L,U,x] = crout_solver(A,b)
disp(A)
disp(b)
disp(L)
disp(U)
disp(L*U)
disp(A*x)

printf("doolittle\n")
[L,U,x] = d_solver(A,b)
disp(A)
disp(b)
disp(L)
disp(U)
disp(L*U)
disp(A*x)
// Factorizacion de Cholesky
function [U, ind] = CholeskyV1(A)
    eps = 1.0e-8
    n = size(A,1)
    U = zeros(n,n)
    for k = 1:n
        if k==1 then
                t = A(k,k)
        else 
                t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
        end
    
        if t <= eps
            printf("Matriz no definida positiva.\n")
            ind = 0
            return
        end
        U(k,k)= sqrt(t)
        for j = k+1:n
            if k==1 then 
                        U(k,j) = A(k,j)/U(k,k)
            else 
                        U(k,j) = ( A(k,j) - U(1:k-1,k)' * U(1:k-1,j) )/U(k,k)
            end
        end
    end
    ind = 1
endfunction

// Factorizacion de Cholesky
function [U,ind] = choleskyV2(A)
    // Factorización de Cholesky.
    // Trabaja únicamente con la parte triangular superior.
    //
    // ind = 1  si se obtuvo la factorización de Cholesky.
    //     = 0  si A no es definida positiva
    //
    //******************
    eps = 1.0e-8
    //******************
    
    n = size(A,1)
    U = zeros(n,n)
    
    t = A(1,1)
    if t <= eps then
        printf('Matriz no definida positiva.\n')
        ind = 0
        return
    end
    U(1,1) = sqrt(t)
    for j = 2:n
        U(1,j) = A(1,j)/U(1,1)
    end
        
    for k = 2:n
        t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
        if t <= eps then
            printf('Matriz no definida positiva.\n')
            ind = 0
            return
        end
        U(k,k) = sqrt(t)
        for j = k+1:n
            U(k,j) = ( A(k,j) - U(1:k-1,k)'*U(1:k-1,j) )/U(k,k)
        end
    end
    ind = 1

endfunction
