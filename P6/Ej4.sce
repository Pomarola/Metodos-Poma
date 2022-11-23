clc
clear

function circ(r,x,y)
    xarc(x-r,y+r,2*r,2*r,0,360*64)
endfunction

//solo circulos
function [radios, centros] = Gers(A)
    [n,m] = size(A)
    centros = diag(A)
    radios = sum(abs(A),'c') - abs(centros)
    
    mx = round (min(centros - radios)-1);
    my = round (min(-radios)-1);
    
    // esquina superior derecha
    
    Mx = round(max(centros+radios)+1);
    My = round(max(radios)+1);
    
    rectangulo = [mx my Mx My];
    
    // dibujamos los autovalores
    plot2d(real(spec(A)),imag(spec(A)),-1,"031","",rectangulo)
    replot(rectangulo);
    xgrid()
    for i=1:n
        circ(radios(i),centros(i), 0)
    end
endfunction

// circulos y autovalores
function [radios, centros] = CircGersValor(A)
    [n,m] = size(A)
    centros = diag(A)
    radios = sum(abs(A),'c') - abs(centros)
    
    mx = round (min(centros - radios)) - 1
    my = round (min(-radios)) - 1
    
    Mx = round (max(centros + radios)) + 1
    My = round (max(radios)) + 1
    
    rect = [mx,my,Mx,My]
    eigen = spec(A)
    plot2d(real(eigen),imag(eigen),-3,"031","",rect)
    xgrid(2,1,7)
    for i=1:n
        circ(radios(i),centros(i), 0)
    end
endfunction

A = [4 1 0; 1 4 1; 0 1 4]
Gers(A)
CircGersValor(A)
