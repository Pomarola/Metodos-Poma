
A = [1 2 -2 1; 4 5 -7 6; 5 25 -15 -3; 6 -12 -6 22]
disp("Arranca Propio")
[L, U, P] = gaussPPLU(A)
disp(A)
disp(L)
disp(U)
disp(P)
disp(P*A)
disp(L*U)
disp("Termina Propio")

// Ejemplos de aplicación

//A = [1 2 -2 1; 4 5 -7 6; 5 25 -15 -3; 6 -12 -6 22]
disp("Ej a")
b = [1 2 0 1]'
//x=
//   9.8333333
//  -6.1666667
//  -5.5
//  -7.5
[x,a] = gaussElimPP(A,b)
disp(x)
disp(a)

disp("Ej b")
b = [2 2 1 0]'

// x = [-1 2 0.00000 1]'
[x,a] = gaussElimPP(A,b)
disp(x)
disp(a)
