function [ NS, DS, J, K ] = Func_N_Square_Nat( xnode , kx, ky, t)
    % NS son las Funciones de Forma simbolicas
    % DS son las derivadas de las Funciones de Forma 
    % N son las Funciones de Forma evaluadas en P
    % D son las derivadas de las Funciones de Forma evaluadas en P
    
    syms e n;
    N1 = (1-e)*(1-n)/4;
    N2 = (1+e)*(1-n)/4;
    N3 = (1+e)*(1+n)/4;
    N4 = (1-e)*(1+n)/4;
    DN1e = diff(N1, e); DN1n = diff(N1, n);
    DN2e = diff(N2, e); DN2n = diff(N2, n);
    DN3e = diff(N3, e); DN3n = diff(N3, n);
    DN4e = diff(N4, e); DN4n = diff(N4, n);

    NS = [N1; N2; N3; N4];
    DS = [DN1e DN1n; DN2e DN2n; DN3e DN3n; DN4e DN4n];    
    
    %% Calculo de la Matriz de Difusion K pero en Coordenadas Naturales
    % int(B' k B t dOmega) = int int B*' k B* t det(J) dn de


    X = simplify(N1 * xnode(1,1) + N2 * xnode(2, 1) + N3 * xnode(3, 1) + N4 * xnode(4, 1));
    Y = simplify(N1 * xnode(1,2) + N2 * xnode(2, 2) + N3 * xnode(3, 2) + N4 * xnode(4, 2));

    DXe = diff(X, e); DXn = diff(X, n);
    DYe = diff(Y, e); DYn = diff(Y, n);

    J = [DXe DYe; DXn DYn];
    detJ = det(J);
    InvJ = inv(J);
    
    % Calculo B* le llamo BM
    DN = [diff(N1, e) diff(N2, e) diff(N3, e) diff(N4, e);
          diff(N1, n) diff(N2, n) diff(N3, n) diff(N4, n)];
    BM = J \ DN;  % inv(J) * DN
    
    k = [kx 0; 0 ky];
    
    M = simplify(BM.' * k * BM * detJ* t);
    
    K = double(int( int(M, e, -1, 1), n, -1, 1));
    
end

