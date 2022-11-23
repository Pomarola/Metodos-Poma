clc
clear

function y = ej8(x)
    y = exp(x) == 3*x
endfunction

function xk1 = g1(xk)
    xk1 = exp(xk) / 3
endfunction

function xk1 = g2(xk)
    xk1 = (exp(xk) - xk) / 2
endfunction

function xk1 = g3(xk)
    xk1 = log(3*xk)
endfunction

function xk1 = g4(xk)
    xk1 = exp(xk) - 2 * xk
endfunction

function r = fijo(ec, g, xk)
    disp(xk)
    if ec(xk) then
        r = xk
    else
        r = fijo(ec, g, g(xk))
    end
endfunction
