function [H,fh] = mat_robin(xi,yi,i,xj,yj,j,h,phi_inf,q,t)
    % Función para calcular la matriz H y el vector f para una arista de un triángulo con condición
    % de contorno Robin.
    % i,j pueden valer 1,2,3 dependiendo el lado... por ejemplo lado 1-2 -> i = 1; j = 2;
    % lado 3-1 -> i = 3; j = 1;
    
    Pi = [xi,yi]; % -> es el punto (i)
    Pj = [xj,yj]; % -> es el punto (j)
    
    lij = norm(Pi-Pj,2);
    fh = zeros(3,1);
    fh([i;j]) = [1;1];
    fh = fh * ((h*phi_inf - q)*lij*t/2);
    
    H = zeros(3,3);
    H([i,j],[i,j]) = [2 1; 1 2];
    H = H * (h*lij*t/6);


end

