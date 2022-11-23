clear
clc

function y = a3(x)
    y = sin(x) - x^2 / 2;
endfunction

function y = b3(x)
    y = exp(-x) - x^4;
endfunction

function y = c3(x)
    y = log10(x) - x + 1;
endfunction

function r = biseccion(f, a, b, tol, iters)
    
    for i = 1:iters
        c = (a + b)/2;
        if (b - c) <= tol then
            break;
        else
            if f(b) * f(c) > 0 then
                b = c;
            else
                a = c;
            end
        end
    end
    if (b - c) > tol then 
        disp("Se alcanzo el maximo de iteraciones");
    end
    r = c;
endfunction

//c3 0.13 0.99
//b3 -1.43 0.81
//a3 0.00 1.39
