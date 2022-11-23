// funcion f dada como string con variable x
// v es el valor a evaluar
// n es el orden de derivacion
// h es el paso de derivacion

function r = derivar(f, v, n, h)
    deff("y=DER0(x)","y="+f);
    if n == 0 then
        r = DER0(v);
    else 
        for i = 1:(n-1)
            deff ("y=DER"+string(i)+"(x)","y=(DER"+string(i-1)+"(x+"+string(h)+")-DER"+string(i-1)+"(x))/"+string(h));
        end
        deff ("y=DERn(x)","y=(DER"+string(n-1)+"(x+"+string(h)+")-DER"+string(n-1)+"(x))/"+string(h));
        r = DERn(v);
    end
endfunction

// usando numderivative con orden 4
function r = derivarNum(f, v, n, h)
    deff("y=DER0(x)","y="+f);
    if n == 0 then
        r = DER0(v);
    else 
        for i = 1:(n-1)
            deff ("y=DER"+string(i)+"(x)","y=numderivative(DER"+string(i-1)+",x,"+string(h)+",4)");
        end
        deff ("y=DERn(x)","y=numderivative(DER"+string(n-1)+",x,"+string(h)+",4)");
        r = DERn(v);
    end
endfunction
