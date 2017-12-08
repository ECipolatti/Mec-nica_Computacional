function [Ke, fe, Kh, fh] = matriz_elemento(elem,xnode,icone,Q,kx,ky,h,phi_inf,qh,t)

    i = icone(elem,1);
    j = icone(elem,2);
    k = icone(elem,3);
    
    xi = xnode(i,1); xj = xnode(j,1); xk = xnode(k,1);
    yi = xnode(i,2); yj = xnode(j,2); yk = xnode(k,2);
    
    [Ke,fe] = mat_diffusion(xi,yi,xj,yj,xk,yk,kx,ky,Q,t);
    [Kh,fh] = mat_robin(xi,yi,i,xj,yj,j,h,phi_inf,qh,t);

end

