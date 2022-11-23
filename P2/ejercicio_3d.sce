function y = hornerDer(p,x)
    n = degree(p);
    a = coeff(p,n);
    b = a;
    for i = 1:n
        if i>1 then b = a + x*b; end
        a = a*x + coeff(p,(n-i));
    end
    y(1) = a
    y(2) = b
endfunction

p = poly([1 2 3],"x","coeff");
x = 3
res = hornerDer(p,x);
printf("horner:                 %e\n", res(1));
printf("derivada por horner:    %e\n", res(2));
