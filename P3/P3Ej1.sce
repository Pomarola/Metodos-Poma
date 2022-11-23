clear
clc

function y = zero(x)
    y = x-x;
endfunction

function y = f1(x)
    y = cos(x) .* cosh(x) + 1
endfunction

x = linspace(-8,8,2000);
plot(x, f1(x), 'r', 'thickness', 2);
plot(x, zero(x), 'b');

disp(f1(1.87));
disp(f1(4.694));
disp(f1(7.8548));
