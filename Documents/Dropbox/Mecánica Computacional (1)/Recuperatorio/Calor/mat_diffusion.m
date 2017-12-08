function [K, fv] = mat_diffusion(xi, yi, xj, yj, xk, yk, kx, ky, Q, t)
    bi = yj - yk;   ci = xk - xj;
    bj = yk - yi;   cj = xi - xk;
    bk = yi - yj;   ck = xj - xi;
    
    A = det([1 xi yi; 1 xj yj; 1 xk yk])/2;
    
    k = [kx 0; 0 ky];
    B = [bi bj bk; ci cj ck]/(2*A);
    
    K = B.' * k * B * t * A;
    
    fv = (Q * t * A/3) * ones(3, 1);

end

