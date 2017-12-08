function [Ni, Nj, Nk, derivadas, var] = ntriang(i,j,k,p,symb)
    xi = i(1); yi = i(2);   % extrae las coordenadas [xi,yi] del punto (i)
    xj = j(1); yj = j(2);   % extrae las coordenadas [xj,yj] del punto (j)
    xk = k(1); yk = k(2);   % extrae las coordenadas [xk,yk] del punto (k)
    xp = p(1); yp = p(2);   % extrae las coordenadas [xp,yp] del punto (p)

    A = det([1 xi yi; 1 xj yj; 1 xk yk])/2; % área del elemento
    % coeficientes de las funciones de forma
    ai = (xj*yk - xk*yj); bi = (yj - yk); ci = (xk - xj);
    aj = (xk*yi - xi*yk); bj = (yk - yi); cj = (xi - xk);
    ak = (xi*yj - xj*yi); bk = (yi - yj); ck = (xj - xi);
    
    % Valores de las Constantes abc por fila 
    var = [ai bi ci; aj bj cj; ak bk ck];
    
    % si se pide simbólico sólo se devuelve las funciones de forma Ni(x,y)
    if symb == 1
        syms x y;
        
        Ni = (ai + bi*x + ci*y)/(2*A);
        Nj = (aj + bj*x + cj*y)/(2*A);
        Nk = (ak + bk*x + ck*y)/(2*A);
        
        %Dirivadas de las funciones de forma
        dxNi = diff(Ni,x);
        dyNi = diff(Ni,y);
        
        dxNj = diff(Nj,x);
        dyNj = diff(Nj,y);
        
        dxNk = diff(Nk,x);
        dyNk = diff(Nk,y);
        derivadas = [dxNi dyNi; dxNj dyNj; dxNk dyNk];
    else % sino se devuelve las funciones de forma evaluadas en el punto (p) [xp,yp]
        Ni = (ai + bi*xp + ci*yp)/(2*A);
        Nj = (aj + bj*xp + cj*yp)/(2*A);
        Nk = (ak + bk*xp + ck*yp)/(2*A);
    end
end
