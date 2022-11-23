clear
clc

function b = ej5(x)
    b = 2*x == 2^x
endfunction

function xk1 = g(xk)
    xk1 = 2^(xk-1)
endfunction

function r = fijo(ec, g, xk)
    if ec(xk) then
        r = xk
    else
        r = fijo(ec, g, g(xk))
    end
endfunction

function r = fijo2(ec, g, xk, iters)
    for i = 1:iters
        if ec(xk) then
            break
        else
            xk = g(xk)
        end
    end
    if ~(ec(xk))
        disp("Se llego al maximo de iteraciones")
    end
    r = xk
endfunction
