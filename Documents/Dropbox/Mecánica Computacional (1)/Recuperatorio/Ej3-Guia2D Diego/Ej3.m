xnode = [0.0 0.0;
         5.0 0.0;
         0.0 5.0;
         5.0 5.0];
 icone = [1 2 3;
          2 4 3];
 % El lado inferior de la placa esta Aislado: q = 0
 % El lado derecho esta a 100 °C
 % El lado superior tiene condicion Robin: h(T - Tref) = 1.2(T - 30)
 % El lado izquierdo tiene flujo: q = 2
 k = 2;         % k del material
 Q = 1.2;       % Carga Q distribuida en el Elemnto 1
 Qp = 5.0;      % Carga Puntual Qp en P = (1, 1) 
 xp = 1;        % Coordenada x de P
 yp = 1;        % Coordenada y de P
 h=1.2;         % h de robin
 Tref = 30;     % T ambiente
 q = 2;         % Flujo en el lado izquierdo
 l = 5;         % placa cuadrada de 5 cm de Lado
 
 K = zeros (size(xnode,1));
 F = zeros (size(xnode,1),1);
 
 % Elemento 1
 ele = icone(1,:);
 coord_ele = xnode(ele,:);
 A = area_ele(coord_ele);
 [Klocal1, Flocal1] = local1(coord_ele,k,k,A,l,q,Qp,xp,yp);
 % Ensamble
 K(ele,ele) = K(ele,ele) + Klocal1;
 F(ele) = F(ele) + Flocal1;
 
 % Elemento 2
 ele = icone(2,:);
 coord_ele = xnode(ele,:);
 A = area_ele(coord_ele);
 [Klocal2, Flocal2] = local2(coord_ele,k,k,A,l,h,Tref,Q);
 % Ensamble
 K(ele,ele) = K(ele,ele) + Klocal2;
 F(ele) = F(ele) + Flocal2;
 
 % Matriz que produce la condicion mixta
 C = (h*l/6)*[0 0 0 0;0 0 0 0;0 0 2 1;0 0 1 2];
 K = K + C;
 
 % Condiciones Dirichlet
 K(2,:)=0; K(2,2)=1; F(2)=100;
 K(4,:)=0; K(4,4)=1; F(4)=100;
 T = K\F;
 
% VISUALIZACION
patch('Faces',icone, 'Vertices',xnode, 'FaceVertexCData',T, ...
    'FaceColor','interp', 'EdgeColor','k', 'EraseMode','normal'); 
% dibujar una marca en los nodos
lin=line(xnode(:,1),xnode(:,2),'LineStyle','none','Marker','o'); 
colorbar;