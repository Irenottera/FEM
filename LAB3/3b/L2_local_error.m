function [ eL2_local ] = L2_local_error (e, ex_sol, Uold, nodes, elements, totnodes, t)

%nodes of the e-th element
cur_nodes = nodes(elements(e, :)',:);

ureal_atnode = zeros(4,1);
unum_atnode = zeros(4,1);

a = size(nodes);
totnodes = a(1);

u_num = Uold(elements(e, :)'); %vector of the numerical solution values at the mesh nodes

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
K = zeros(8, 2);
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
    K(i*2-1:i*2, :) = [Y(i), 0 ; 0, X(i)];
end

%compute the function in the 4 integration points and then sum them.

w1 = w1*J(1);
w2 = w2*J(2);
w3 = w3*J(3);
w4 = w4*J(4);


for(i=1:4)
    ureal_atnode(i) = ex_sol(X(i), Y(i), t);
end

for (i=1:4) %PROBABLY THE MISTAKE IS HERE
    for (j=1:4)
    unum_atnode(i) = unum_atnode(i)+ u_num(j)*PHI(i,j)  ; 
    end
end


        
        
        k1 = (abs(ureal_atnode(1)-unum_atnode(1)))^2;
        k2 = (abs(ureal_atnode(2)-unum_atnode(2)))^2;
        k3 = (abs(ureal_atnode(3)-unum_atnode(3)))^2;
        k4 = (abs(ureal_atnode(4)-unum_atnode(4)))^2;
        eL2_local = w1*k1+w2*k2+w3*k3+w4*k4 ;
        eL2_local = sqrt(eL2_local);
        



end










