function y = hornerR(p,x,n)
    y = coeff(p,n);
    if n<>degree(p) then
        y = y + x*hornerR(p,x,n+1);
    end
endfunction

function y = hornerRec(p,x)
    y = hornerR(p,x,0);
endfunction

function y = hornerIte(p,x)
    n = degree(p);
    y = coeff(p,n);
    for i = 0:(n-1)
        y = y*x + coeff(p,(n-1-i));
    end
endfunction

p = poly([1 2 3],"x","coeff");
printf("%e\n", hornerRec(p,3));
printf("%e\n", horner(p,3));
printf("%e", hornerIte(p,3));
