function [K] = mat_stiff_triang_nat(xi,yi,xj,yj,xk,yk,E,nu,t)
    syms e n;

    N1 = 1 - e - n;
    N2 = e;
    N3 = n;

    x = N1*xi + N2*xj + N3*xk;
    y = N1*yi + N2*yj + N3*yk;

    D = mat_stress(E,nu);   % OJO: fijarse si es tensión o deformación

    I = [1 0 0 0;
        0 0 0 1;
        0 1 1 0];

    DN1 = [diff(N1,e), 0, diff(N2,e), 0, diff(N3,e), 0]; % u con respecto a ksi
    DN2 = [diff(N1,n), 0, diff(N2,n), 0, diff(N3,n), 0]; % u con respecto a eta
    DN3 = [0, diff(N1,e), 0, diff(N2,e), 0, diff(N3,e)]; % v con respecto a ksi
    DN4 = [0, diff(N1,n), 0, diff(N2,n), 0, diff(N3,n)]; % v con respecto a eta

    DN = [DN1; DN2; DN3; DN4];

    J = [diff(x,e) diff(y,e);
        diff(x,n) diff(y,n)];

    S = [inv(J), zeros(2,2);
        zeros(2,2) inv(J)];

    detJ = det(J);

    B = I * S * DN;
    M = B.' * D * B * t * detJ;

    K = double(int(int(M,e,0,1-n),n,0,1));
   
end
