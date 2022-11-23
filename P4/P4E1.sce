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
