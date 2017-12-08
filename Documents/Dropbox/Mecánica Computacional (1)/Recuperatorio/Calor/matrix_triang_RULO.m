function [K, var, A] = matrix_triang_RULO(Pi,Pj,Pk,kx,ky)
    xi = Pi(1); xj = Pj(1); xk = Pk(1);
    yi = Pi(2); yj = Pj(2); yk = Pk(2);

    ai = xj*yk - xk*yj;
    aj = xk*yi - xi*yk;
    ak = xi*yj - xj*yi;
        
    bi = yj - yk;
    bj = yk - yi;
    bk = yi - yj;

    ci = xk - xj;
    cj = xi - xk;
    ck = xj - xi;
    
    var = [ai aj ak; bi bj bk; ci cj ck];

    D = [1 xi yi; 1 xj yj; 1 xk yk];

    A = det(D)/2;

    K = [bi*bi*kx + ci*ci*ky, bi*bj*kx + ci*cj*ky, bi*bk*kx + ci*ck*ky; ...
        bj*bi*kx + cj*ci*ky, bj*bj*kx + cj*cj*ky, bj*bk*kx + cj*ck*ky; ...
        bk*bi*kx + ck*ci*ky, bk*bj*kx + ck*cj*ky, bk*bk*kx + ck*ck*ky]/(4*A);
end

