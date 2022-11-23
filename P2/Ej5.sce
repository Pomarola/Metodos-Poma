function r = taylor(f, v, n)
    h = 0.01;
    deff("y=DF0(x)","y="+f);
    coeffs(1) = DF0(0)
    for i = 1:n
        coeffs(i+1) = derivarNum(f, 0, i, h) / factorial(i);
    end
    p = poly(coeffs, "x", "coeff");
    disp(p);
    r = horner(p, v);
endfunction

//for c=1:3
//    printf("Exponencial: %e\n", exp(1));
//    printf("Taylor de exp: %e\n", taylors("exp(x)", 1, c));
//end
