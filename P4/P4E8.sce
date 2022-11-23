// Ejemplos de aplicación

A = [1.012 -2.132 3.104; -2.132 4.096 -7.013; 3.104 -7.013 0.014]
disp("Arranca Propio")
[L, U, P] = gaussPPLU(A)
disp(A)
disp(L)
disp(U)
disp(P)
disp(P*A)
disp(L*U)
disp("Termina Propio")

disp("Arranca Scilab")
[L,U,P] = lu(A)
disp(A)
disp(L)
disp(U)
disp(P)
disp(P*A)
disp(L*U)
disp("Termina Scilab")

A = [2.1756 4.0231 -2.1732 5.1967; -4.0231 6.0000 0 1.1973; -1.0000 5.2107 1.1111 0 ; 6.0235 7.0000 0 4.1561]

disp("Arranca Propio")
[L, U, P] = gaussPPLU(A)
disp(A)
disp(L)
disp(U)
disp(P)
disp(P*A)
disp(L*U)
disp("Termina Propio")

disp("Arranca Scilab")
[L,U,P] = lu(A)
disp(A)
disp(L)
disp(U)
disp(P)
disp(P*A)
disp(L*U)
disp("Termina Scilab")

//Dan lo mismo, ejecutar P4E7 antes de esto
