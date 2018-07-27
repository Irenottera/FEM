function [ Ke ] = LocalstiffK( e, nodes, elements, k)

%memorize the nodes of the e-th element
cur_nodes = nodes(elements(e, :)',:);

a = size(nodes);
totnodes = a(1);

%integration points 
int_points = zeros (4 , 2);
int_points(1,1)= -1/sqrt(3);
int_points(1,2)= -1/sqrt(3);
int_points(2,1)= 1/sqrt(3);
int_points(2,2)= -1/sqrt(3);
int_points(3,1)= 1/sqrt(3);
int_points(3,2)= 1/sqrt(3);
int_points(4,1)= -1/sqrt(3);
int_points(4,2)= 1/sqrt(3);

w1 = 1;
w2 = 1;
w3 = 1;
w4 = 1;


%basis functions at integration points 
PHI = zeros(4,4);
%e-derivative of basis functions at integration points 
DPHIDE = zeros(4,4);
%n-derivative of basis functions at integration points 
DPHIDN = zeros(4,4);
%x(e,n) and y(e,n) evaluated at integration points 
X = zeros(4,1);
Y = zeros(4,1);
DXDE = zeros(4,1);
DXDN = zeros(4,1);
DYDE = zeros(4,1);
DYDN = zeros(4,1);
J = zeros(4,1);
DEDX = zeros(4, 1);
DEDY = zeros(4, 1);
DNDX = zeros(4, 1);
DNDY = zeros(4, 1);
DPHIDX = zeros(4,4);
DPHIDY = zeros(4,4);
B = zeros(4,1);
F = zeros(4,1);
Fe = zeros(totnodes,1);
K = zeros(4, 1);
Ke = zeros(totnodes,totnodes);

for i=1:4
    for j=1:4
    PHI(i,j) = master_shape(j, int_points(i,1), int_points(i,2));
    DPHIDE(i,j) = master_shape_deriv_x(j, int_points(i,1), int_points(i,2));
    DPHIDN(i,j) = master_shape_deriv_y(j, int_points(i,1), int_points(i,2));
    X(i) = X(i) + cur_nodes(j,1)*master_shape(j,int_points(i,1),int_points(i,2));
    Y(i) = Y(i) + cur_nodes(j,2)*master_shape(j,int_points(i,1),int_points(i,2));
    DXDE(i) = DXDE(i) + cur_nodes(j,1)*master_shape_deriv_x(j,int_points(i,1),int_points(i,2));
    DXDN(i) = DXDN(i) + cur_nodes(j,1)*master_shape_deriv_y(j,int_points(i,1),int_points(i,2)); 
    DYDE(i) = DYDE(i) + cur_nodes(j,2)*master_shape_deriv_x(j,int_points(i,1),int_points(i,2));
    DYDN(i) = DYDN(i) + cur_nodes(j,2)*master_shape_deriv_y(j,int_points(i,1),int_points(i,2)); 
    end
end

for (i=1:4)
    J(i) = DXDE(i)*DYDN(i)-DXDN(i)*DYDE(i);
    DEDX(i) = (1/J(i))*DYDN(i);
    DEDY(i) = -(1/J(i))*DXDN(i);
    DNDX(i) = -(1/J(i))*DYDE(i);
    DNDY(i) = (1/J(i))*DXDE(i);
end

for (i=1:4)
    for (j=1:4)
    DPHIDX(i,j) = DPHIDE(i,j)*DEDX(i)+DPHIDN(i,j)*DNDX(i);
    DPHIDY(i,j) = DPHIDE(i,j)*DEDY(i)+DPHIDN(i,j)*DNDY(i);
    end
end

for (i=1:4)
    K(i) = k(X(i), Y(i));
    %K(i*2-1:i*2, :) = [Y(i), 0 ; 0, X(i)];
end

%compute the function in the 4 integration points and then sum them. do the
%same fot the load vector. 

w1 = w1*J(1);
w2 = w2*J(2);
w3 = w3*J(3);
w4 = w4*J(4);

ii=1;
jj=1;

for (i=elements(e, :))

    jj = 1;

    for (j=elements(e, :))        
        
        v1 = [DPHIDX(1,ii), DPHIDY(1,ii)]';
        v2 = [DPHIDX(2,ii), DPHIDY(2,ii)]';
        v3 = [DPHIDX(3,ii), DPHIDY(3,ii)]';
        v4 = [DPHIDX(4,ii), DPHIDY(4,ii)]';
        vv1 = [DPHIDX(1,jj), DPHIDY(1,jj)]';
        vv2 = [DPHIDX(2,jj), DPHIDY(2,jj)]';
        vv3 = [DPHIDX(3,jj), DPHIDY(3,jj)]';
        vv4 = [DPHIDX(4,jj), DPHIDY(4,jj)]';
       
        k1 = dot(K(1).*v1, vv1);
        k2 = dot(K(2).*v2, vv2);
        k3 = dot(K(3).*v3, vv3);
        k4 = dot(K(4).*v4, vv4);
        Ke(i,j) = w1*k1+w2*k2+w3*k3+w4*k4 ;
        
        jj = jj+1;  
    end
ii = ii+1;
end



end

