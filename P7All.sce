clc
clear

// 2 puntos lineal, 3 puntos cuadratica, 4 puntos cubica
// x = poly(0,"x")
// para acotar el error la derivada enesima de f, aplicado a lo que daria el maximo valor en el ej 1 es 0.6

// Polinomio interpolador de Lagrange
// Parametros: punto variable p, lista de x0, lista de y0
function w = Lagrange(p,x,y)
        w = 0
        n = length(x)
    for i=1:n do
        w = w + L(p,i,x)*y(i)
    end
endfunction

// Función L_i(x) del polinomio interpolador de Lagrange
// Parametros: punto variable p, indice i, lista de y0
function w = L(p,i,x)
    w = 1
    n = length(x)
    for j=1:n do
        if j<>i then
            w = w*(p-x(j))/(x(i)-x(j))
        end
    end
endfunction

// Diferencias Divididas de newton
// Parametros: punto variable p, lista de x0, lista de y0
function w = DD_Newton(p,x,y)
    w = 0
    n = length(x)
    for j=n:-1:2
        w = (w + DD(x(1:j),y(1:j)))*(p-x(j-1))
    end
    w = w + y(1)
endfunction


// Diferencias divididas
function w = DD(x,y)
    n = length(x)
    if n==2 then
        w = (y(n)-y(1))/(x(n)-x(1))
    else
        w = (DD(x(2:n),y(2:n))-DD(x(1:n-1),y(1:n-1)))/(x(n)-x(1))
    end
endfunction

// Método de Diferencias Divididas de Newton
// Formula con multiplicaciones encajadas
// NO FUNCIONA ARREGLAR
function w = DD_NewtonME(x,y)
    // Entrada: x,y = vectores puntos de interpolación (x,y)
    // Salida: w = polinomio de diferencias divididas de Newton
    s = poly(0,"x")
    n = length(x)
    w =  DD(x,y)
    for j=n-1:-1:1
        w = DD(x(1:j),y(1:j)) + (s-x(j))*w
    end
endfunction

// Error de interpolación
function w = err(p,x,cot)
    // Entrada: p = valor real, x = nodos de interpolación, cot = cota de |f^(n))|
    // Salida: w = error de interpolación en x = p
    n = length(x)
    w = cot/(factorial(n))
    for i=1:n do
        w = w*abs(p - x(i))
    end
endfunction







function [x,a] = gaussElimPP(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana con pivoteo.  

[nA,mA] = size(A) 
[nb,mb] = size(b)

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;

a = [A b]; // Matriz aumentada

// Eliminación progresiva
n = nA;
for i = 1:n-1
    pivot = i
    amax = abs(a(i,i))  //pivoteo
    for j=i+1:n
        if abs(a(j,i)) > amax then
            pivot = j
            amax = a(j,i)
        end
    end
    temp = a(pivot,:)
    a(pivot,:) = a(i,:)
    a(i,:) = temp
    
    for j=i+1:n
        mji = a(j,i)/a(i,i)
        a(j,i) = 0
        a(j,i+1:n+1) = a(j,i+1:n+1) - mji * a(i,i+1:n+1)
    end
end

x(n) = a(n,n+1)/a(n,n);
for i = n-1:-1:1
    x(i) = ( a(i,n+1) - ( a(i, i+1:n) * x(i+1:n) ) ) / a(i,i)
end

endfunction

// Minimos cuadrados

// Aproximación polinomial de mínimos cuadrados polinomial para matrices con rango completo
function [p,err] = MinCuad_pol(A,b)
    // Entrada: b = vectores 1xn
    // Salida: p = polinomio de mínimos cuadrados; err = vector de errores (eps = Ax-b)
     [w,a] = gaussElimPP((A')*A,(A')*(b'))
     p = poly(w,"x","coeff")
     err = A*w-b'
endfunction

// Matriz del método de mínimo cuadrados polinomial
function A = matriz_mc(x,n)
    // Entrada: x,y = vectores 1xn; n = grado de la aproximación
    // Salida: A = matriz del método de mínimo cuadrados
    m = length(x)
    A = ones(m,1)
    for j=2:(n+1) do
        A = [A,(x').^(j-1)]
    end
endfunction

// Funcion que calcula los nodos de raices segun el grado.
function r = Cheb(n)
    for k=0:n-1
        r(k+1) = cos(%pi/2*(1+2*k)/n)
    end
endfunction

// Recibe un numero n
// Devuelve el polinomio de chebyshev de ese grado con sus raíces
function [T,r] = Chebyshev(n)
    t(1) = 1
    t(2) = poly([0],"x","r")
    for i = 3:n+1
        t(i) = poly([0 2], "x", "coeff")*t(i-1)-t(i-2) 
    end
    T = t(n+1)
    r = roots(T)
endfunction

// Raices de Chebyshev en cualquier intervalo
function x = NodosChebyshev(n,a,b)
    [pol,r] = Chebyshev(n)
    for i = 1 : n
        x(i) = ((b+a) + r(i) * (b - a))/2
    end
endfunction




// EJERCICIO 1 ----------------------------------------------
/*
x1 = [0.2,0.4]
y1 = [1.2214,1.4918]

xc = [0,0.2,0.4,0.6]
yc = [1.0,1.2214,1.4918,1.8221]

p = poly(0,"x")

// Calculando polinomios
lagrangeLineal = Lagrange(p,x1,y1)
lagrangeCubico = Lagrange(p,xc,yc)

disp("Resultado de evaluar polinomio lienal en 1/3")
disp(horner(lagrangeLineal,1/3))
// = 1.4016667
disp("Cota de error lineal:")
disp(abs((1/3-0.2)*(1/3-0.4)*1/2*%e^(0.6)))
// = 0.0080983
disp("Error del polinomio lineal")
disp(abs(1.395612425-horner(lagrangeLineal,1/3)))
// = 0.0060542

disp("Resultado de evaluar polinomio cubico en 1/3")
disp(horner(lagrangeCubico,1/3))
// = 1.3955494
disp("Cota de error cubico:")
disp((1/3-0)*(1/3-0.2)*(1/3-0.4)*(1/3-0.6)*1/factorial(4)*%e^0.6)
// = 0.00006
disp("Error del polinomio cubico")
disp(abs(1.395612425-horner(lagrangeCubico,1/3)))
// = 0.000063

x1 = [0.2,0.4]
y1 = [1.2214,1.4918]

xc = [0,0.2,0.4,0.6]
yc = [1.0,1.2214,1.4918,1.8221]

p = poly(0,"x")

newtonLineal = DD_Newton(p,x1,y1)
newtonCubica = DD_Newton(p,xc,yc)

disp("Newton")

disp("Resultado de evaluar polinomio lienal en 1/3")
disp(horner(newtonLineal,1/3))
// = 1.4016667
disp("Cota de error lineal:")
disp(abs((1/3-0.2)*(1/3-0.4)*1/2*%e^(0.6)))
// = 0.0080983
disp("Error del polinomio lineal")
disp(abs(1.395612425-horner(newtonLineal,1/3)))
// = 0.0060542

disp("Resultado de evaluar polinomio cubico en 1/3")
disp(horner(newtonCubica,1/3))
// = 1.3955494
disp("Cota de error cubico:")
disp((1/3-0)*(1/3-0.2)*(1/3-0.4)*(1/3-0.6)*1/factorial(4)*%e^0.6)
// = 0.00006
disp("Error del polinomio cubico")
disp(abs(1.395612425-horner(newtonCubica,1/3)))
// = 0.000063
*/
// EJERCICIO 2 ----------------------------------------------

//f(x) polinomio de orden menor o igual a n, demostrar que el polinomio de lagrange de orden n para f es exacto
//
//vemos la forma del error del polinomio de interpolacion de orden n
//f(x) - p(x) = (x-x0)(x-x1)...(x-xn) / (n+1)!  f^(n+1) (eps(x))
//
//Como f es un polinomio de grado a lo sumo n, la derivada n+1-esima de f es 0. Entonces tenemos f(x) - p(x) = 0   => f(x) = p(x) 

// EJERCICIO 3 ----------------------------------------------

// EJERCICIO 4 ----------------------------------------------
/*
x = [2,2.1,2.2,2.3,2.4,2.5]
y = [0.2239,0.1666,0.1104,0.0555,0.0025,-0.0484]

p = poly(0,"x")

p = DD_Newton(p,x,y)
disp("El polinomio interpolante es: ")
disp(p)

w1 = horner(p,2.15)
err1 = err(2.15,x,1) // |j_0'(x)| <= 1
disp("El valor aproximado de J_0(2.15) es: "+string(w1))
disp("con error: "+string(err1)+" < 0.5D-06")

w2 = horner(p,2.35)
err2 = err(2.35,x,1) // |j_0'(x)| <= 1
disp("El valor aproximado de J_0(2.35) es: "+string(w2))
disp("con error: "+string(err2)+" < 0.5D-06")
*/
// EJERCICIO 5 ----------------------------------------------
/*
// Despejamos los valores para x=0,1,2
// Utilizando los polinomios de interpolacion que son dato
// f(0) = 1, f(1) = 3 y f(2) = 3

// Despejamos f(3)
x = [1 2 3]

p = poly(0,"x")

L1 = L(2.5,1,x)
L2 = L(2.5,2,x)
L3 = L(2.5,3,x)

c1 = L1*3
c2 = L2*3
c3 = L3

a = (3-c1-c2)/c3

x = [0 1 2 3]
y = [1 3 3 a]

res = Lagrange(2.5,x,y)
disp(res)

*/
// EJERCICIO 6 ----------------------------------------------
/*
// DATOS
DDF1 = 2 // f[-1]
DDF2 = 1 // f[-1,1]
DDF3 = -2 // f[-1,1,2]
DDF4 = 2 // f[-1,1,2,4]

x = [-1 1 2 4]

// y0 = DDF1
// DDF2 = (y1 - y0) / (x1 - x0) => DDF2*(x1-x0) + y0 = y1

y1 = DDF2*(x(2)-x(1))+DDF1

// DDF3 = f[x1,x2]-f[x0,x1] / (x2-x0)
// f[x1,x2] = y2 - y1 / x2 - x1
// DDF3*(x2-x0) + DDF2 = (y2 - y1) / (x2 - x1)
// [DDF3*(x2-x0)+DDF2]*(x2-x1) + y1

y2 = (DDF3*(x(3)-x(1)) + DDF2)*(x(3)-x(2)) + y1

function y = pol(x)
    y = DDF1 + (x+1)*(DDF2 + (x-1)*(DDF3 + (x-2)*DDF4))
endfunction
printf("Aproximacion de f(0): %.2f\n", pol(0))

err0 = err(0,x,33.6)

printf("Cota de error de aproximacion de f(0): %f\n", err0)
// El resultado es muy extremo, revisar
*/

// EJERCICIO 7 ----------------------------------------------
/*
x = [0 0.15 0.31 0.5 0.6 0.75]
y = [1 1.004 1.31 1.117 1.223 1.422]

// Grado 1
A = matriz_mc(x,1)
[p1,err1] = MinCuad_pol(A,y)
disp(p1)
A = matriz_mc(x,2)
[p2,err2] = MinCuad_pol(A,y)
disp(p2)
A = matriz_mc(x,3)
[p3,err3] = MinCuad_pol(A,y)
disp(p3)

rango = [-0:0.0001:1]
plot(rango, horner(p1, rango),"r")
plot(rango, horner(p2, rango),"g")
plot(rango, horner(p3, rango),"b")
plot(x', y', "r*")
h1 = legend(["Grado 1", "Grado 2", "Grado 3"])
*/

// EJERCICIO 8 ----------------------------------------------

/*
x = [4 4.2 4.5 4.7 5.1 5.5 5.9 6.3 6.8 7.1]
y = [102.56 113.18 130.11 142.05 167.53 195.14 224.87 256.73 299.5 326.72]

disp("ítem a)")

disp("(#) n=1.")
A = matriz_mc(x,1)
deter = det((A')*A)
disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
disp("El polinomio de mínimos cuadrados de grado 1 es:")
[p1,err1] = MinCuad_pol(A,y)
disp(p1)

disp("(#) n=2.")
A = matriz_mc(x,2)
deter = det((A')*A)
disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
disp("El polinomio de mínimos cuadrados de grado 2 es:")
[p2,err2] = MinCuad_pol(A,y)
disp(p2)

disp("(#) n=3.")
A = matriz_mc(x,3)
deter = det((A')*A)
disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
disp("El polinomio de mínimos cuadrados de grado 3 es:")
[p3,err3] = MinCuad_pol(A,y)
disp(p3)


disp("(#) Analizamos los errores err = norm(Ax-y,2).")
disp("Para la aproximación lineal: "+string(norm(err1,2)))
disp("Para la aproximación cuadrática: "+string(norm(err2,2)))
disp("Para la aproximación cúbica: "+string(norm(err3,2)))
disp("Podemos decir que es mejor la aproximación cúbica en este caso.")

inf = min(x)-1
sup = max(x)+1
rango = [inf:0.001:sup]
plot(rango, horner(p1, rango),"r")
plot(rango, horner(p2, rango),"g")
plot(rango, horner(p3, rango),"b")
plot(x', y', "r*")
h1 = legend(["Grado 1", "Grado 2", "Grado 3"])
*/

// EJERCICIO 9 ----------------------------------------------
/*
function y=f(x)
    y = 1./(1+x.^2)
endfunction

sup = 5
inf = -sup
rango = [inf:0.0001:sup]
func = f(rango)

for n=2:2:14
    if (n <> 8 && n <> 12)
        nodos = linspace(inf,sup,n+1)
        p = poly(0,"x")
        p = Lagrange(p,nodos,f(nodos))
        aprox = horner(p,rango)
        plot(rango,aprox,"b")
        plot(rango, func,"r")
        plot(rango, abs(aprox-func),"g")
        h1 = legend(["Aprox Grado "+ string(n), "1/1+x^2", "Error"])
        figure
        sleep(1,"s")
    end
end
*/

// EJERCICIO 10 ----------------------------------------------
/*
rango = [-1:0.0001:1]
func = exp(rango)

nodos = Cheb(4)
p = poly(0,"x")
p = DD_Newton(p,nodos',exp(nodos'))


plot(rango,horner(p,rango),"r")
plot(rango,func,"b")
plot(rango,abs(horner(p,rango)-func),"g")
h1 = legend(["Aprox Grado 3", "exp(x)", "Error"])
a=get("current_axes")//get the handle of the newly created axes
a.axes_visible="on"; // makes the axes visible
a.y_location= "middle"
*/

// EJERCICIO 11 ----------------------------------------------
/*
a=0
b=%pi/2
rango = [a:0.00001:b]
func = cos(rango)
grado = 3

nodos = NodosChebyshev(grado+1,a,b)

p = poly(0,"x")
p = DD_Newton(p,nodos',cos(nodos'))

plot(rango,horner(p,rango),"r")
plot(rango,func,"b")
plot(rango,abs(horner(p,rango)-func),"g")
h1 = legend(["Aprox Grado 3", "cos(x)", "Error"])
*/
