close all; clear all;

%% DATOS

E = 2e6;
nu = 0.2;
t = 0.5;

icone = [ 1 2 3 4;
          2 5 6 3 ];

xnode = [ 0 0;
          1 0;
          1 1;
          0 1;
          2 0;
          2 1 ];
      
      
%% ENSAMBLE
M = size(icone,1);  % cantidad de elementos
N = size(xnode,1);	% cantidad de nodos

nrowk = [];
ncolk = [];
coefk = [];

for e = 1:M
    %% CÁLCULO DE LA MATRIZ ELEMENTAL
    i = icone(e,1);
    j = icone(e,2);
    k = icone(e,3);
    l = icone(e,4);
    
    xi = xnode(i,1); xj = xnode(j,1); xk = xnode(k,1); xl = xnode(l,1); 
    yi = xnode(i,2); yj = xnode(j,2); yk = xnode(k,2); yl = xnode(l,2);
    
    Ke = mat_stiff_square_nat(xi,yi,xj,yj,xk,yk,xl,yl,E,nu,t);
    
    %% ENSAMBLE
    for r = 1 : 4
        for c = 1 : 4
            ncolk = [ncolk, 2*icone(e,c)-1];  %#ok<*AGROW>
            ncolk = [ncolk, 2*icone(e,c)];
            ncolk = [ncolk, 2*icone(e,c)-1];
            ncolk = [ncolk, 2*icone(e,c)];
            
            
            nrowk = [nrowk, 2*icone(e,r)-1];
            nrowk = [nrowk, 2*icone(e,r)-1];
            nrowk = [nrowk, 2*icone(e,r)];
            nrowk = [nrowk, 2*icone(e,r)];
            
            coefk = [coefk, Ke(2*r-1,2*c-1)];
            coefk = [coefk, Ke(2*r-1,2*c)];
            coefk = [coefk, Ke(2*r,2*c-1)];
            coefk = [coefk, Ke(2*r,2*c)];
        end
    end
end

K = sparse(nrowk,ncolk,coefk);
K = full(K);