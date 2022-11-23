clear
clc

function [L, U] = dooLittle(A)
    n = size(A, 1)
    U(1,:) = A(1,:)
    L = eye(n,n)
    L(2:n,1) = A(2:n,1)/U(1,1)
    U(2, 2:n) = A(2,2:3) - (L(2,2) * U(1,2:3))
    disp(U)
    disp(L)
    for i = 2:n
        for j = i+1:n
            disp(L(j,1:i-1))
            disp(U(1:i-1,i))
            L(j,i) = (A(j,i) - (L(j,1:i-1) * U(1:i-1,i))) / U(i,i)
            U(j,i) = A(j,i) - L(j,1:j-1) * U(1:j-1,i)
        end
    end
endfunction

A = [1 2 3 4; 1 4 9 16; 1 8 27 64; 1 16 81 256]
[L, U] = dooLittle(A)
disp(A)
disp(L)
disp(U)
disp(L*U)

A = [1 2 -2 1; 4 5 -7 6; 5 25 -15 -3; 6 -12 -6 22]
[L, U] = dooLittle(A)
disp(A)
disp(L)
disp(U)
disp(L*U)
