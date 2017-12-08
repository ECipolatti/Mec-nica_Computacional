function [K, B, D, A] = mat_stiff_triang_cart(xi,yi,xj,yj,xk,yk,E,nu,t)
   % matriz de rigidez para un elemento triangular
   A = det([1 xi yi; 1 xj yj; 1 xk yk])/2;
   
   bi = yj - yk; ci = xk - xj;
   bj = yk - yi; cj = xi - xk;
   bk = yi - yj; ck = xj - xi;

   B = [bi 0 bj 0 bk 0; 
        0 ci 0 cj 0 ck; 
        ci bi cj bj ck bk]/(2*A);
    
   D = mat_stress(E,nu);   % OJO: fijarse si es tensión o deformación
   
   K = B.' * D * B * t * A;
   
end

