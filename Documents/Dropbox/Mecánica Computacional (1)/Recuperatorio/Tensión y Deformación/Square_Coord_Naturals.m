xnode = [0.0 0.0
         2.2 0.5
         2.0 1.0
         0.0 1.0 ];
P = [0.5 0.5];
kx = 1;
ky = 1;
t = 0.1;
[ NS, DS, J, K ] = Func_N_Square_Nat( xnode, kx, ky, t )