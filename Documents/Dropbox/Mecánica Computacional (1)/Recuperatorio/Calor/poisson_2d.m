close all; clear all;

kx = 1; ky = 1;
Q = 0;

t = 1;
% Matriz de conexión (indica qué nodos forman cada elemento)

icone = [1, 2, 4;
         2, 5, 4;
         2, 3, 5;
         4, 5, 6];
% Posición (x,y) de cada nodo de la malla
% xnod = [0 0;2.5 0; 2.5 2.5; 5 0; 5 2.5; 5 5];
xnod = [0 0;
        0.5 0;
        1 0;
        0.25, sin(pi/3)/2;
        0.75, sin(pi/3)/2;
        0.5, sin(pi/3)];

% 0- Dirichlet; 1- Neumann; 2- Robin
% B: [elemento,tipo,valor,h,phi_inf]
% e.g. Dirichlet: [4,0,0,0,0] -> elemento: 4, tipo: Dirichlet, valor: 0[, h: 0, phi_inf: 0]
% e.g. Neumann:   [3,1,2,0,0] -> elemento: 3, tipo: Neumann,   valor: 2[, h: 0, phi_inf: 0]
% e.g. Robin:     [5,2,0,2,5] -> elemento: 5, tipo: Robin,     valor: 0,  h: 2, phi_inf: 5
% B = [4,0,0; 5,0,0; 6,0,0];
% B = [1,0,100; 2,0,100; 3,0,100; 4,0,100; 6,0,100; 7,0,500; 8,0,500; 9,0,500];
DIR = [1 100; 2 100; 3 100];
NEU = [1 4 100; 4 6 100];
ROB = [3 5 0 10 50; 5 6 0 10 50]; 
Puntual = [];  % Elemento , x, y, carga

% [phi,K,f,KB,fB] = fem_2d(kx,ky,Q,icone,xnod,B);
[phi,K,f,KB,fB] = fem_2d(kx,ky,Q,Puntual,t,icone,xnod,DIR,NEU,ROB)

