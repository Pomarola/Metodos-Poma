clear
clc

function y = f3(x)
    y = x ^ 2 / 4 - sin(x);
endfunction

function r = secante(f, x0, x1, tol, iters)
    for i = 1:iters
        aux = x1 - f(x1) * ((x1 - x0) / (f(x1) - f(x0)));
        x0 = x1;
        x1 = aux;
        if abs(x1 - x0) < tol then break;
        end
    end
    if abs(x1 - x0) > tol then
        disp("Se alcanzo el maximo de iteraciones");
    end
    r = x1;
endfunction

// f3 1.93
// f3 0
