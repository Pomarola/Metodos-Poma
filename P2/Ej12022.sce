clc
clear

function r = robusta(p)
    c = coeff(p,0);
    b = coeff(p,1);
    a = coeff(p,2);
    d = sqrt(b^2-4*a*c);
    if b < 0 then
        r(1) = (2*c) / (-b + d) 
        r(2) = (-b + d)/(2*a);
    else
        r(1) = (-b - d)/(2*a);
        r(2) = (2*c) / (-b - d)
    end
endfunction

p = poly([-0.0001 10000.0 0.0001],"x","coeff");
e1 = 1e-8;
roots1 = robusta(p);
r1 = roots1(1);
r2 = roots1(2);
printf("robusto: %e %e \n", r1, r2);
roots2 = roots(p);
printf("scilab:  %e %e\n", roots2(1), roots2(2));
error = abs(r2-e1)/e1;
printf("error de robusto para la segunda raiz : %e)\n", error);
