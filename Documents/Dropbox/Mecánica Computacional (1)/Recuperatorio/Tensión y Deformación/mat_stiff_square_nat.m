function [K, J] = mat_stiff_square_nat(xi,yi,xj,yj,xk,yk,xl,yl,E,nu,t)
    syms e n;

    Ni = (1-e)*(1-n)/4;
    Nj = (1+e)*(1-n)/4;
    Nk = (1+e)*(1+n)/4;
    Nl = (1-e)*(1+n)/4;

    x = Ni*xi + Nj*xj + Nk*xk + Nl*xl;
    y = Ni*yi + Nj*yj + Nk*yk + Nl*yl;

    D = mat_stress(E,nu);       % OJO: fijarse si es tensión o deformación

    I = [1 0 0 0;
        0 0 0 1;
        0 1 1 0];

    DN1 = [diff(Ni,e), 0, diff(Nj,e), 0, diff(Nk,e), 0, diff(Nl,e), 0]; % u con respecto a ksi
    DN2 = [diff(Ni,n), 0, diff(Nj,n), 0, diff(Nk,n), 0, diff(Nl,n), 0]; % u con respecto a eta
    DN3 = [0, diff(Ni,e), 0, diff(Nj,e), 0, diff(Nk,e), 0, diff(Nl,e)]; % v con respecto a ksi
    DN4 = [0, diff(Ni,n), 0, diff(Nj,n), 0, diff(Nk,n), 0, diff(Nl,n)]; % v con respecto a eta

    DN = [DN1; DN2; DN3; DN4];

    J = [diff(x,e) diff(y,e);
        diff(x,n) diff(y,n)];

    S = [inv(J), zeros(2,2);
        zeros(2,2) inv(J)];

    detJ = det(J);

    B = I * S * DN;
    M = simplify(B.' * D * B * t * detJ);

    K = double(int(int(M, e, -1, 1), n, -1, 1 ));
   
end
