function x = sistCholesky(A, b)
    [U, ind] = cholesky(A)
    g = progresive_sust(U', b)
    x = regresive_sust(U, g)
endfunction

D = [16 -12 8; -12 18 -6; 8 -6 8]
b = [76 -66 46]'
x = sistCholesky(D, b)
disp(x)


