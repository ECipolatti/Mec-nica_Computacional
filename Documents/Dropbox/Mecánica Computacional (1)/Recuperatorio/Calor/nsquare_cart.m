function [ Ni, Nj, Nk, Nl, derivadas ] = nsquare_cart(i, j , k, l, p, symb)
    xi = i(1); yi = i(2);   % extrae las coordenadas [xi,yi] del punto (i)
    xj = j(1); yj = j(2);   % extrae las coordenadas [xj,yj] del punto (j)
    xk = k(1); yk = k(2);   % extrae las coordenadas [xk,yk] del punto (k)
    xl = l(1); yl = l(2);   % extrae las coordenadas [xl,yl] del punto (l)
    xp = p(1); yp = p(2);   % extrae las coordenadas [xp,yp] del punto (p)

    a = (xj - xi)/2;
    b = (yk - yj)/2;
    
    xp = xp - a;
    yp = yp - b;
    
    % si se pide simbólico sólo se devuelve las funciones de forma Ni(x,y)
    if symb == 1
        syms x y;
        
        Ni = (a - x)*(b - y)/(4*a*b);
        Nj = (a + x)*(b - y)/(4*a*b);
        Nk = (a + x)*(b + y)/(4*a*b);
        Nl = (a - x)*(b + y)/(4*a*b);
        
        % Dirivadas de las funciones de forma
        dxNi = diff(Ni,x);  dyNi = diff(Ni,y);
        dxNj = diff(Nj,x);  dyNj = diff(Nj,y);
        dxNk = diff(Nk,x);  dyNk = diff(Nk,y);
        dxNl = diff(Nl,x);  dyNl = diff(Nl,y);
        
        derivadas = [dxNi dyNi; dxNj dyNj; dxNk dyNk; dxNl dyNl];
        
    else % sino se devuelve las funciones de forma evaluadas en el punto (p) [xp,yp]
        syms x y;
        Ni = (a - x)*(b - y)/(4*a*b);
        Nj = (a + x)*(b - y)/(4*a*b);
        Nk = (a + x)*(b + y)/(4*a*b);
        Nl = (a - x)*(b + y)/(4*a*b);
        
        dxNi = diff(Ni,x);  dyNi = diff(Ni,y);
        dxNj = diff(Nj,x);  dyNj = diff(Nj,y);
        dxNk = diff(Nk,x);  dyNk = diff(Nk,y);
        dxNl = diff(Nl,x);  dyNl = diff(Nl,y);
        
        Ni = (a - xp)*(b - yp)/(4*a*b);
        Nj = (a + xp)*(b - yp)/(4*a*b);
        Nk = (a + xp)*(b + yp)/(4*a*b);
        Nl = (a - xp)*(b + yp)/(4*a*b);
        
        derivadas = [subs(dxNi, y, yp) subs(dyNi, x, xp);
                     subs(dxNj, y, yp) subs(dyNj, x, xp);
                     subs(dxNk, y, yp) subs(dyNk, x, xp);
                     subs(dxNl, y, yp) subs(dyNl, x, xp)];
    end
end

