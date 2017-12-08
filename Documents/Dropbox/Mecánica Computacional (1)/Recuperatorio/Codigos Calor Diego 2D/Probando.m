xnode = [0.0 0.0; 5.0 0.0; 5.0 5.0; 0.0 5.0];
kx = 2;
ky = 2;

%% Matriz Difusiva K para Triangulos y Cuadrangulos en Coor. Cartesianas
[K1] = MatCondCar(xnode,kx,ky)

%% Matriz Difusiva K para Triangulos y Cuadrangulos en Coor. Naturales(simb)
[K2] = MatCondNat(xnode,kx,ky)

%% Matriz Difusiva K para Triangulos y Cuadrangulos en Coor. Naturales (Integracion Numérica)
[K3] = MatCondNatIN(xnode,kx,ky)