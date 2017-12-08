close all; clear all;

a = 1;
b = 1;
t = 0.50;

x1 = 0; x2 = a; x3 = a; x4 = 0;
y1 = 0; y2 = 0; y3 = b; y4 = b;

E = 2e6;
nu = 0.2;

syms e n;

N1 = (1-e)*(1-n)/4;
N2 = (1+e)*(1-n)/4;
N3 = (1+e)*(1+n)/4;
N4 = (1-e)*(1+n)/4;

x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4;

D = mat_stress(E,nu);

I = [1 0 0 0;
     0 0 0 1;
     0 1 1 0];

DN1 = [diff(N1,e), 0, diff(N2,e), 0, diff(N3,e), 0, diff(N4,e), 0]; % u con respecto a ksi 
DN2 = [diff(N1,n), 0, diff(N2,n), 0, diff(N3,n), 0, diff(N4,n), 0]; % u con respecto a eta 
DN3 = [0, diff(N1,e), 0, diff(N2,e), 0, diff(N3,e), 0, diff(N4,e)]; % v con respecto a ksi 
DN4 = [0, diff(N1,n), 0, diff(N2,n), 0, diff(N3,n), 0, diff(N4,n)]; % v con respecto a eta 

DN = [DN1; DN2; DN3; DN4];

J = [diff(x,e) diff(y,e);
     diff(x,n) diff(y,n)];
 
S = [inv(J), zeros(2,2);
     zeros(2,2) inv(J)];
 
detJ = det(J);
 
B = I * S * DN;
M = B.' * D * B * t * detJ;

K = double(int(int(M,e,-1,1),n,-1,1));
disp(K);