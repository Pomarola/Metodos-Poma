clear
clc

function r = sist(v)
    x = v(1)
    y = v(2)
    fx = 1 + x^2 - y^2 + exp(x) * cos(y)
    fy = 2*x*y + exp(x) * sin(y)
    r = [fx; fy]
endfunction

function r = newton(f, xn, n)
    if n == 0 then
        r = xn
    else
        j = numderivative(f,xn)^-1
        xn1 = xn - j * f(xn)
        r = newton(f, xn1, n-1)
    end
endfunction

printf("%f\n", newton(sist, [-1; 4], 5))

function r = newton2(f, v0, n)
    for i = 1:n
        j = numderivative(f, v0)^-1
        v1 = v0 - j * f(v0)
        v0 = v1
    end
    r = v1
endfunction

printf("%f\n", newton2(sist, [-1; 4], 5))
