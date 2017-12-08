function [Ni, Nj, Nk, derivadas, var] = ntriang_elem(elem,icone,xnode,xp,yp,symb)
    i = icone(elem,1);     % nodo (i) del elemento
    j = icone(elem,2);     % nodo (j) del elemento
    k = icone(elem,3);     % nodo (k) del elemento
    
    % coordenadas [x,y] de cada nodo
    xi = xnode(i,1); xj = xnode(j,1); xk = xnode(k,1);
    yi = xnode(i,2); yj = xnode(j,2); yk = xnode(k,2);

    A = det([1 xi yi; 1 xj yj; 1 xk yk])/2; % área del elemento
    
    % coeficientes de las funciones de forma
    ai = (xj*yk - xk*yj); bi = (yj - yk); ci = (xk - xj);
    aj = (xk*yi - xi*yk); bj = (yk - yi); cj = (xi - xk);
    ak = (xi*yj - xj*yi); bk = (yi - yj); ck = (xj - xi);
    
    var = [ai aj ak; bi bj bk; ci cj ck];
    
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
        
        derivadas = [];
    end
end
