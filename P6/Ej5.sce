clc
clear

//potencia(A, z0, max_iter). Dada una matriz A y un vector z0 (estimación de posible autovector) y una cantidad de iteraciones max_iter
// aproxima al autovalor cuyo módulo es el radio espectral (mayor valor absoluto).
function rho = potencia(A, z0, max_iter)
    sz = size(z0, 1)
    for i = 1:max_iter
        w = A*z0
        //disp(w)    
        if(i <> max_iter)
            z = w / norm(w, %inf)
            z0 = z
        end
    end
    //Elegimos la componente de mayor valor absoluto
    k = 1
    for i = 2:sz
        if( abs(w(i)) > abs(w(k)))
            k = i
        end
    end
    rho = w(k) / z0(k)
endfunction

function [v, zn, iters]= mpotencia(A,z0,eps,maxiter)
    V = max(abs(spec(A)))
    v = 0
    iters = 1
    w = A*z0
    zn = w/norm(w)
    [m,j] = max(abs(w))
    v = m/z0(j) // Falta checkear que z0(j) != 0
    er1 = abs(V-v)
    er2 = norm(zn-z0)
    err = max(er1,er2)
    while (iters <= maxiter && err > eps)
        z0 = zn
        w = A*z0
        zn = w/norm(w,'inf')
        [m,j] = max(abs(w))
        v = m/z0(j) // Falta checkear que z0(j) != 0
        er1 = abs(V-v)
        er2 = norm(zn-z0,'inf')
        err = max(er1,er2)
        iters = iters+1
    end
    printf("iteraciones: %d\n", iters)
endfunction

// ----------------- EJERCICIO 5 -----------------

A1 = [6 4 4 1; 4 6 1 4; 4 1 6 4; 1 4 4 6]
A2 = [12 1 3 4; 1 -3 1 5; 3 1 6 -2; 4 5 -2 -1]
Z0 = ones(4,1)
eps = 10^(-8)
maxiter = 1000
[v1,z1, iters1] = mpotencia(A1,Z0,eps,maxiter)
printf("Matriz A1\n")
printf("Autovalor dominante: %f\n", v1)
printf("Autovector asociado:")
disp(z1)
printf("\n")

[v2,z2, iters2] = mpotencia(A2,Z0,eps,maxiter)
printf("Matriz A2\n")
printf("Autovalor dominante: %f\n", v2)
printf("Autovector asociado:")
disp(z2)
