function [phi,K,f,KB,fB] = fem_2d(kx,ky,Q,P,t,icone,xnode,DIR,NEU,ROB)
    % kx: conductividad en sentido x
    % ky: conductividad en sentido y
    % Q: fuente
    % t: espesor
    % icone: matriz de conectividad de elementos triangulares
    % xnod: posición (x,y) de cada nodo de la malla
    % DIR: matriz con los datos de los nodos con condición Dirichlet
    % NEU: matriz con los datos de los bordes con condición Neumann
    % ROB: matriz con los datos de los bordes con condición Robin
    
    M = size(icone,1);  % cantidad de elementos
    N = size(xnode,1);	% cantidad de nodos

    nrowk = [];
    ncolk = [];
    coefk = [];
    
    nrowf = [];
    coeff = [];
    
    for e = 1:M
        i = icone(e,1);
        j = icone(e,2);
        k = icone(e,3);

        xi = xnode(i,1); xj = xnode(j,1); xk = xnode(k,1);
        yi = xnode(i,2); yj = xnode(j,2); yk = xnode(k,2);

        [Ke,fe] = mat_diffusion(xi,yi,xj,yj,xk,yk,kx,ky,Q,t);
        
        for r = 1 : 3
            for c = 1 : 3
                ncolk = [ncolk, icone(e,c)];
                nrowk = [nrowk, icone(e,r)]; %#ok<*AGROW>
                coefk = [coefk, Ke(r,c)];
            end
            nrowf = [nrowf, icone(e,r)];
            coeff = [coeff, fe(r)];
        end
        
    end
    
    %% Cálculo de cargas puntuales
    NP = size(P,1);  
    for n = 1 : NP
        e = P(n,1);         % elemento al que se aplica la carga puntual
        i = icone(e,1);     % nodo (i) del elemento
        j = icone(e,2);     % nodo (j) del elemento
        k = icone(e,3);     % nodo (k) del elemento

        % coordenadas [x,y] de cada nodo
        xi = xnode(i,1); xj = xnode(j,1); xk = xnode(k,1);
        yi = xnode(i,2); yj = xnode(j,2); yk = xnode(k,2);
        
        xp = P(n,2);    % coordenada x del punto donde se aplica la carga
        yp = P(n,3);    % coordenada y del punto donde se aplica la carga
        G = P(n,4);     % valor de la carga aplicada
        
        % funciones de forma evaluadas en el punto [xp,yp]
        [Ni,Nj,Nk] = ntriang([xi,yi],[xj,yj],[xk,yk],[xp,yp],0);
        % cálculo de la carga distribuida para este elemento. luego hay que restarla ya que
        % previamente se calculó para todos los elementos una carga distribuida
        [~,fe] = mat_diffusion(xi,yi,xj,yj,xk,yk,kx,ky,Q,t);
        
        % cálculo de la carga puntual
        fp = G*t*[Ni;Nj;Nk];
       
        % se actualizan los vectores de coeficientes para el armado del vector f (sparse)
        for r = 1 : 3
            nrowf = [nrowf, icone(e,r)];
            coeff = [coeff, fp(r)-fe(r)]; %ELIMINAR Fe si queremos la carga distribuida
                                          %sobre el elemento donde est� la carga puntual;
        end
    end
    
    %% se ensamblan la matriz y el vector como sparse
    K = sparse(nrowk,ncolk,coefk);
    f = sparse(nrowf,1,coeff);

    
    %% Elementos de los Bordes
    
    % NEUMANN
    NB = size(NEU,1);
    for n = 1 : NB
        Pi = xnode(NEU(n,1),:);
        Pj = xnode(NEU(n,2),:);
        lij = norm(Pi-Pj,2);
        q = NEU(n,3);
        
        f(NEU(n,1)) = f(NEU(n,1)) - q*lij*t/2;
        f(NEU(n,2)) = f(NEU(n,2)) - q*lij*t/2;
    end
    
    %% ROBIN
    NB = size(ROB,1);
    for n = 1 : NB
        Pi = xnode(ROB(n,1),:);
        Pj = xnode(ROB(n,2),:);
        lij = norm(Pi-Pj,2);
        q = ROB(n,3);
        h = ROB(n,4);
        phi_inf = ROB(n,5);
        K(ROB(n,1),ROB(n,1)) = K(ROB(n,1),ROB(n,1)) + 2*(h*lij*t/6);
        K(ROB(n,1),ROB(n,2)) = K(ROB(n,1),ROB(n,2)) + (h*lij*t/6);
        K(ROB(n,2),ROB(n,1)) = K(ROB(n,2),ROB(n,1)) + (h*lij*t/6);
        K(ROB(n,2),ROB(n,2)) = K(ROB(n,2),ROB(n,2)) + 2*(h*lij*t/6);
        
        f(ROB(n,1)) = f(ROB(n,1)) + (h*phi_inf - q)*lij*t/2;
        f(ROB(n,2)) = f(ROB(n,2)) + (h*phi_inf - q)*lij*t/2;
    end
    
    
    KB = K;
    fB = f;
        
    %% DIRICHLET
    NB = size(DIR,1);     % cantidad de elementos de borde Dirichlet
    for n = 1 : NB
        idx = DIR(n,1);
        val = DIR(n,2);
        fB(idx) = val; %#ok<*SPRIX>
        KB(idx,:) = zrow(idx,N);
    end
    
    %% Resolución del sistema
    phi = KB\fB;        % resuelve el sistema
    KB = full(KB);      % pasa de sparse a full la matriz alterada para condiciones de borde
    K = full(K);        % pasa de sparse a full la matriz original (real para el problema)
    f = full(f);        % pasa de sparse a full el vector RHS (real para el problema)
    fB = full(fB);      % pasa de sparse a full el vector RHS alterado para condiciones de borde
    phi = full(phi);    % pasa de sparse a full el vector solución
end

