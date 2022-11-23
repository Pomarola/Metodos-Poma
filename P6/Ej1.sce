clc
clear

function [c,r] = cotasGers(A)
    n = size(A,1)
    c = diag(A)
    r = sum(abs(A),'c') - abs(c)
    
    for i = 1:n
        mprintf("Circulo centrado en %f de radio %f \n", c(i), r(i))
    end
endfunction

A = [1 0 0; -1 0 1; -1 -1 2]
// circulo centrado en (1,0) radio 0
// circulo centrado en (0,0) radio 2
// circulo centrado en (2,0) radio 2
// los autovalores estan en el intervalo [-2,4]
cotasS = spec(A)
disp("Autovalores A:")
disp(cotasS)
cotasGers(A)
 
B = [1 0 0; -0.1 0 0.1; -0.1 -0.1 2]
// circulo centrado en (1,0) radio 0
// circulo centrado en (0,0) radio 0.2
// circulo centrado en (2,0) radio 0.2
// los autovalores estan en el intervalo [-0.2,0.2] U {1} U [1.8,2.2]
cotasS = spec(B)
disp("Autovalores B:")
disp(cotasS)
cotasGers(B)

C = [1 0 0; -0.25 0 0.25; -0.25 -0.25 2]
// circulo centrado en (1,0) radio 0
// circulo centrado en (0,0) radio 0.5
// circulo centrado en (2,0) radio 0.5
// los autovalores estan en el intervalo [-0.5,0.5] U {1} U [1.5,2.5]
cotasS = spec(C)
disp("Autovalores C:")
disp(cotasS)
cotasGers(C)

D = [4 -1 0; -1 4 -1; -1 -1 4]
// circulo centrado en (4,0) radio 1
// circulo centrado en (4,0) radio 2
// circulo centrado en (4,0) radio 2
// los autovalores estan en el intervalo [2,6]
cotasS = spec(D)
disp("Autovalores D:")
disp(cotasS)
cotasGers(D)

E = [3 2 1; 2 3 0; 1 0 3]
// circulo centrado en (3,0) radio 3
// circulo centrado en (3,0) radio 2
// circulo centrado en (3,0) radio 1
// los autovalores estan en el intervalo [0,6]
cotasS = spec(E)
disp("Autovalores E:")
disp(cotasS)
cotasGers(E)

F = [4.75 2.25 -0.25; 2.25 4.75 1.25; -0.25 1.25 4.75]
// circulo centrado en (4.75,0) radio 2.5
// circulo centrado en (4.75,0) radio 3.5
// circulo centrado en (4.75,0) radio 1.5
// los autovalores estan en el intervalo [1.25,8.25]
cotasS = spec(F)
disp("Autovalores F:")
disp(cotasS)
cotasGers(F)


