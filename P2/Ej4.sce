// funcion f es la ley de la función dada por un string, usa como 
// variable x
// v es el valor donde se evaluará la derivada
// n es el orden de derivación
// h es el paso de derivación

function valor = derivada(f,v,n,h)
    deff("y=DF0(x)","y="+f);
    if n==0 then valor = DF0(v);
    else
        for i=1:(n-1)
        deff("y=DF"+string(i)+"(x)","y=(DF"+string(i-1)+"(x+"+string(h)+")-DF"+string(i-1)+"(x))/"+string(h));
        end
        deff("y=DFn(x)","y=(DF"+string(n-1)+"(x+"+string(h)+")-DF"+string(n-1)+"(x))/"+string(h));
        valor = DFn(v);
    end
endfunction

// usando numderivative
// esta función utiliza un orden para numderivative igual a 4
function valor = derivadaNum(f,v,n,h)
    deff("y=DF0(x)","y="+f);
    if n==0 then valor = DF0(v);
    else
        for i=1:(n-1)
        deff("y=DF"+string(i)+"(x)","y=numderivative(DF"+string(i-1)+",x,"+string(h)+",4)");
        end
        deff("y=DFn(x)","y=numderivative(DF"+string(n-1)+",x,"+string(h)+",1)");
        valor = DFn(v);
    end
endfunction


function value = TaylorSwift (f, n, v)
    deff("y=DF0(x)","y="+f);
    t(1) = DF0(0)
    for i=1:n
        t(i+1) = derivadaNum(f,0,i,0.01)/factorial(i);
    end
    p = poly(t,"x","coeff");
    disp(p);
    value = horner(p,v)
    
endfunction


for i=1:8
    printf("Exponencial: %e\n", exp(1));
    printf("Taylor de exp: %e\n", TaylorSwift("exp(x)", i, 1));
end


