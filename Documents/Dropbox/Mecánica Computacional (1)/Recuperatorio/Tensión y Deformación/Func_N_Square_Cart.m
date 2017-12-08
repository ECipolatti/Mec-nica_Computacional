function [ NS, DS, N, D ] = Func_N_Square_Cart( xnode, P )
    % NS son las Funciones de Forma simbolicas
    % DS son las derivadas de las Funciones de Forma 
    % N son las Funciones de Forma evaluadas en P
    % D son las derivadas de las Funciones de Forma evaluadas en P
    % Función de Forma Genérica para Cuadrangulos
    % Ni = ai + bi x + ci y + di xy
    
    syms x y
    A = [1, xnode(1,1), xnode(1,2), xnode(1,1)*xnode(1,2);
         1, xnode(2,1), xnode(2,2), xnode(2,1)*xnode(2,2);
         1, xnode(3,1), xnode(3,2), xnode(3,1)*xnode(3,2);
         1, xnode(4,1), xnode(4,2), xnode(4,1)*xnode(4,2)];
    Alpha = zeros(4, 4);
    for i=1 : 4
        f = zeros(4,1);
        f(i) = 1;
        Alpha(i, :) = (A\f)';
    end
    
    N1 = Alpha(1, 1) + Alpha(1, 2) * x + Alpha(1, 3) * y + Alpha(1, 4) * (x * y);
    N2 = Alpha(2, 1) + Alpha(2, 2) * x + Alpha(2, 3) * y + Alpha(2, 4) * (x * y);
    N3 = Alpha(3, 1) + Alpha(3, 2) * x + Alpha(3, 3) * y + Alpha(3, 4) * (x * y);
    N4 = Alpha(4, 1) + Alpha(4, 2) * x + Alpha(4, 3) * y + Alpha(4, 4) * (x * y);
    
    DN1x = diff(N1, x); DN1y = diff(N1, y);
    DN2x = diff(N2, x); DN2y = diff(N2, y);
    DN3x = diff(N3, x); DN3y = diff(N3, y);
    DN4x = diff(N4, x); DN4y = diff(N4, y);
    
    NS = [N1; N2; N3; N4]; 
    DS = [DN1x DN1y; DN2x DN2y; DN3x DN3y; DN4x DN4y];
    
    N = [subs(N1, [x,y], P); 
             subs(N2, [x, y], P);
             subs(N3, [x, y], P);
             subs(N4, [x, y], P)];
    
    D = [subs(DN1x, y, P(2)) subs(DN1y, x, P(1));
            subs(DN2x, [x, y], P) subs(DN2y, [x, y], P);
            subs(DN3x, [x, y], P) subs(DN3y, [x, y], P); 
            subs(DN4x, [x, y], P) subs(DN4y, [x, y], P)];        

end

