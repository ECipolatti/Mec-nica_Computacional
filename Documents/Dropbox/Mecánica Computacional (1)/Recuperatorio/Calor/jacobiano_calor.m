close all; clear all;

a = 0.40;
b = 0.40;
t = 0.10;

x1 = 0; x2 = a; x3 = a; x4 = 0;
y1 = 0; y2 = 0; y3 = b; y4 = b;

kx = 1;
ky = 1;

E = 2*10^6;
nu = 0.2;

syms e n x y;

N1 = (1-e)*(1-n)/4;
N2 = (1+e)*(1-n)/4;
N3 = (1+e)*(1+n)/4;
N4 = (1-e)*(1+n)/4;

x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4;

DNe = [ diff(N1,e), diff(N2,e), diff(N3,e), diff(N4,e)];
DNn = [ diff(N1,n), diff(N2,n), diff(N3,n), diff(N4,n)];
DN = [DNe; DNn];

J = [diff(x,e) diff(y,e); diff(x,n) diff(y,n)];
detJ = det(J);

D = mat_stress(E,nu);

B1 = J\DN;             % inv(J)*DN
k = [kx 0; 0 ky];
M = simplify(B1.' * k * B1 * detJ);

K = int(int(M,e,-1,1),n,-1,1)