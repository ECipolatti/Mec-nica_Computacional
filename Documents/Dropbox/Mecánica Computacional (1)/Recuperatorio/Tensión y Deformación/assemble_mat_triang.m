close all; clear all;

E = 2.1e5;
nu = 0.3;
t = 10;

icone = [ 1 2 4;
          2 3 4 ];

xnode = [ 0 0;
          1000 0;
          1000 500;
          0 500];

M = size(icone,1);  % cantidad de elementos
N = size(xnode,1);	% cantidad de nodos

nrowk = [];
ncolk = [];
coefk = [];

for e = 1:M
    i = icone(e,1);
    j = icone(e,2);
    k = icone(e,3);
    
    xi = xnode(i,1); xj = xnode(j,1); xk = xnode(k,1);
    yi = xnode(i,2); yj = xnode(j,2); yk = xnode(k,2);
    
    Ke = mat_stiff_triang_nat(xi,yi,xj,yj,xk,yk,E,nu,t);
    
    for r = 1 : 3
        for c = 1 : 3
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