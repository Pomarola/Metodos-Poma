function [L, U, P] = gaussPPLU(A)

n = size(A,1) 
U = A
L = eye(n,n)
P = eye(n,n)
// Eliminación progresiva
for i = 1:n-1
    
    pivot = i
    amax = abs(U(i,i))  //pivoteo
    for j=i+1:n
        if abs(U(j,i)) > amax then
            pivot = j
            amax = U(j,i)
        end
    end
    // Swap U
    temp = U(pivot,i:n)
    U(pivot,i:n) = U(i,i:n)
    U(i,i:n) = temp
    
    // Swap L
    temp = L(pivot,1:i-1)
    L(pivot,1:i-1) = L(i,1:i-1)
    L(i,1:i-1) = temp
    
    // Swap P
    temp = P(pivot,:)
    P(pivot,:) = P(i,:)
    P(i,:) = temp
    
    
    for j= i+1:n
        L(j,i) = U(j,i)/U(i,i)
        U(j,i) = 0
        U(j,i+1:n) = U(j,i+1:n) - L(j,i) * U(i,i+1:n)
    end
end

endfunction


// Ejemplos de aplicación

A = [2 1 1 1; 4 3 3 1; 8 7 9 5; 6 7 9 8]

disp(A)
[L, U, P] = gaussPPLU(A)
disp(L)
disp(U)
disp(P)
disp(P*A)
disp(L*U)
