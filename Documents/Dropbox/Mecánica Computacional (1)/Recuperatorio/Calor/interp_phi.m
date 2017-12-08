function [phi_p] = interp_phi(xp,yp,xi,yi,phi_i,xj,yj,phi_j,xk,yk,phi_k)
    [Ni,Nj,Nk] = ntriang(xi,yi,xj,yj,xk,yk,xp,yp,0);
    
    phi_p = Ni*phi_i + Nj*phi_j + Nk*phi_k;
end

