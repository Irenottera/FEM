% %% CODING PROJECT 2 - Irene Natale - in2582 
% The square domain [0,1]x[0,1] has 4 different boundaries. I numbered them in a anti clockwise order from x=0.
% The nodes are numbered from (0,0) as first node to (1,1) as last node. The elements are numbered in the same way. 
% Also the nodes for each element are numbered in the same way from 1 to 4.
% The problem is defined on the following data 
% q = @(x,y) (x+y)*cos(x-y);
% gd_2 = @(x,y) cos(x);
% gd_4 = @(x,y) cos(x-1);
% gn_1 = @(x,y) 0;
% gn_3 = @(x,y) sin(1-y);
% the pedix of the functions points on which boundary is referred. q is the forcing function. 
% I used square elements and basis functions of the type a_1 + a_2 * x + a_3 * y + a_4 * x * y = 0
% The integrals of the finite element methods have been computed with a Gaussian Quadrature. For 2d integrals
% over each element I used 4 integration points 
% (-1/sqrt(3), -1/sqrt(3))
% (1/sqrt(3), -1/sqrt(3))
% (1/sqrt(3), 1/sqrt(3))
% (-1/sqrt(3), 1/sqrt(3))
% derived from a double 1d integration. Weight functions are 1 in both 1d and 2d integrations.
% 
% The following code is able to solve non homogeneus elliptic partial differential equations in two dimentions with 
% spacially varying Dirichlet BCs on the left and right boundary and spacially varying Neumann BCs on the upper and lower boundary. 

clc 
clear all 

%%%%%%%%%%%%%%DATA DEFINITION
q = @(x,y) (x+y)*cos(x-y);
gd_2 = @(x,y) cos(x);
gd_4 = @(x,y) cos(x-1);
gn_1 = @(x,y) 0;
gn_3 = @(x,y) sin(1-y);
ex_sol = @(x,y) cos(x-y);


%definition of first and last node of the grid 
x1 = 0;
y1 = 0;
xn = 1;
yn = 1;
%number of nodes (per row, per columns, total)
nnodesx = 4;
nnodesy = 4;
totnodes = nnodesx*nnodesy;
%number of elements (per ro, per column, total)
nelemx = nnodesx - 1;
nelemy = nnodesy - 1;
totelem = nelemx*nelemy;
%vector of nodes with their corrispective coordinates
[ nodes ] = RNODE( nnodesx, nnodesy, x1, y1, xn, yn );

%matrix in which elements are linked to their corrispective nodes 
[ elements ] = RELEM( totelem, nelemx );

%creation of the total stiffness matrix K and load vector F
K = zeros (totnodes, totnodes);
F = zeros (totnodes, 1);
for i=1:totelem
    [Ke, Fe] = Localstiff( i, nodes, elements, q );
    K = K + Ke;
    F = F + Fe;
end

%%%%%%%%%%%%%%BOUNDARY CONDITIONS 

%%%%NEUMANN

%Selection of elements I need to put the BCs on
ii=1;
for (i=1:totelem)
    if (mod(i, nelemx) == 0)
        cur_elem(ii, :) = elements(i, :);
        indexelem(ii) = i;
        ii = ii+1;
    end
end

%Gamma vector with the Neuman BCs increments 
G = zeros(totnodes, 1);
for (i=indexelem)    
    Ge = Boundary_Neum(i, nelemx, nelemy, nodes, elements, gn_1, gn_3);    
    G = G + Ge;    
end

%update the load vector 
Fnew = F - G;

%%%%DIRICHLET

%Add Dirichlet BCs
Fdi = zeros(totnodes, 1);
Fdi = Boundary_Dirch_F(Fnew, K, nelemx, nelemy, nodes, elements, gd_2, gd_4);

%Reduction of K and load vector
Kdiri = K  (    nnodesx+1:totnodes-(nnodesx)      ,   nnodesy+1:totnodes-(nnodesy) );
Fdiri = Fdi ( nnodesx+1:totnodes-(nnodesx) );
 
%%%%%%%%%%%%%%SOLUTION ON THE SYSTEM

%Conjugate gradient solver
U = pcg(Kdiri,Fdiri, 1e-12, 1000);
%U = [gd_2( nodes(nnodesx:-1:1 , 1) , 0)   ]

%Solution in nodes-matrix
ii = 1;
U1 = zeros(nnodesx, nnodesy);
for (i = 2: (nnodesy-1) )
    for (j=1:nnodesx)
        U1(i, j) = U(ii);
        ii = ii+1;
    end
end
%Add bounday conditions values
U1 = [(gd_4(nodes(totnodes:-1:totnodes-nnodesx+1 , 1)))' ; U1(2: (nnodesy-1) ,1:nnodesx ) ; (gd_2( nodes(nnodesx:-1:1 , 1) , 0))'];

%%%%%%%%%%%%%%REAL SOLUTION 

%matrix of the real solution values
Ureal = zeros(nnodesx, nnodesy);
ii=1;
vector = [nnodesy:-1:1];
for (i = 1:nnodesy)
    for(j=1:nnodesx)        
        Ureal(i,j) = ex_sol(nodes(ii, 1), nodes(ii,2));
        ii = ii+1;
    end
end

%%%%%%%%%%%%%%PLOTS

figure
subplot(1, 2, 1)
mesh(Ureal)
title('real solution')
subplot(1,2,2)
mesh (U1)
title('numerical solution')
zlim([0.5, 1])

%%%%%%%%%%%%%%ERROR COMPUTATION 

eL2=0;
for (i = 1:nnodesx)
    for (j = 1:nnodesy)
        eL2 = eL2 + ((Ureal(i,j) - U1(i,j)))^2;
    end
end
eL2 = sqrt(eL2);
eL2


%% 

function [ nodes ] = RNODE( nnodesx, nnodesy, x1, y1, xn, yn )
%this function creates a vector with the coordinates of all mesh nodes
nelemx = nnodesx - 1;
nelemy = nnodesy - 1;

dx = (xn-x1)/nelemx;
dy = (yn-y1)/nelemy;
totnodes = nnodesx*nnodesy;

nodes = zeros(totnodes, 2);

for (i=1:totnodes)
    
    if (mod(i, nnodesx) ~= 0 ) 
        nodes(i, 1)= (mod(i, nnodesx) -1)*dx ;
        nodes(i, 2)= (floor(i/nnodesx))*dy;
    elseif  (mod(i, nnodesx) == 0)
        nodes(i, 1)= xn ;
        nodes(i, 2)= nodes(i-1, 2);
    end
    
end

end

%%

function [ elements ] = RELEM( totelem, nelemx )
%this fucntion creates a map of the mesh. The matrix that it creates points
%out which nodes are on which element.

elements = zeros(totelem, 4);
for (i = 1: totelem) %number of the element
        if (mod(i, nelemx)~=0)
            elements(i, 1) =  i+ floor(i/nelemx);
            elements(i, 2) =  elements(i, 1) + 1 ;
            elements(i, 3) =  i+ floor(i/nelemx) + nelemx +2;
            elements(i, 4) =  i+ floor(i/nelemx) + nelemx +1;
        elseif ( mod(i, nelemx)==0)
            elements(i, 1) =  i+ floor(i/nelemx)-1;
            elements(i, 2) =  elements(i, 1)+1 ;
            elements(i, 3) =  i+ floor(i/nelemx) + nelemx +1;
            elements(i, 4) =  i+ floor(i/nelemx) + nelemx ;
        end
    
end

end

%%
%this function creates local stiffnes matrix and load vector for a
%specified element of the mesh
function [ Ke, Fe ] = Localstiff( e, nodes, elements, q )

%memorize the nodes of the e-th element
cur_nodes = zeros (4, 2);
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


%basis functions at itnegration points 
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
    F(i) = q(int_points(i,1), int_points(i,2));
    K(i*2-1:i*2, :) = [X(i), 0 ; 0, Y(i)];
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
       
        k1 = dot(K(1:2, :)*v1, vv1);
        k2 = dot(K(3:4, :)*v2, vv2);
        k3 = dot(K(5:6, :)*v3, vv3);
        k4 = dot(K(7:8, :)*v4, vv4);
        Ke(i,j) = w1*k1+w2*k2+w3*k3+w4*k4 ;
        
        jj = jj+1;  
    end
ii = ii+1;
end

ii=1;

for (i=elements(e, :))
    
      %i is the node
      %ii is the index of the basis function (1:4)
      f1 =  q(X(1), Y(2))*PHI(1, ii);
      f2 =  q(X(2), Y(2))*PHI(2, ii);
      f3 =  q(X(3), Y(3))*PHI(3, ii);
      f4 =  q(X(4), Y(4))*PHI(4, ii);      
      Fe(i) = w1*f1+w2*f2+w3*f3+w4*f4;
      ii = ii+1;
   
end




end

%% 
function [ Ge ] = Boundary_Neum( e, nelemx, nelemy, nodes, elements, gn_1, gn_3)
% e know neumann bcs are defined only on the right boundary of
% the domain, so the element of this vector will be different from zero
% only on the nodes-index that actually are on this boundary. 

% this function constructs the local vector for a specific element e

totelem = nelemx*nelemy;
totnodes = (nelemx+1)*(nelemy+1);
% I find the elements that ly on the right boundary 
cur_elem = zeros(nelemy, 4);   
cur_nodes = zeros (4, 2);
cur_nodes = nodes(elements(e, :)',:);  
  
%2 integration points  - remember it's a 1-d integral
int_points = zeros (2,2);
int_points(1, :) = [1, -1/sqrt(3)];
int_points(2, :) = [1, 1/sqrt(3)];
%weight functions 
w1 = 1;
w2 = 1;

%basis functions on the 2 integration points 
PHI=zeros(2,4);
X = zeros(2,1);
Y = zeros(2,1);
GN = zeros(2,1);
DXDN = zeros(2,1);
DYDN = zeros(2,1);
J = zeros (2,1);
Ge = zeros(totnodes, 1);
for (i = 1:2)
    for (j = 1:4)
        PHI(i,j) = master_shape(j, int_points(i, 1), int_points(i,2));
        X(i) = X(i) + cur_nodes(j,1)*master_shape(j,int_points(i,1),int_points(i,2));
        Y(i) = Y(i) + cur_nodes(j,2)*master_shape(j,int_points(i,1),int_points(i,2));
        DXDN(i) = DXDN(i) + cur_nodes(j,1)*master_shape_deriv_y(j,int_points(i,1),int_points(i,2)); 
        DYDN(i) = DYDN(i) + cur_nodes(j,2)*master_shape_deriv_y(j,int_points(i,1),int_points(i,2));
    end
end

for (i=1:2)
    GN(i) = gn_3(X(i), Y(i));
end

for (i=1:2)
    J(i) = sqrt( ( DXDN(i) )^2 + ( DYDN(i) )^2 );
end

ii=1;
for (i=elements(e, :))
    Ge(i) = GN(1)*PHI(1,ii)*J(1)*w1 + GN(2)*PHI(2,ii)*J(2)*w2;
    ii = ii+1;
end




end

%% 
function [ Fe ] = Boundary_Dirch_F(F, K, nelemx, nelemy, nodes, elements, gd_2, gd_4)
%this function modifies the load vector adding the Dirichlet boudnary
%conditions
 
totnodes = (nelemx+1)*(nelemy+1);
Fe = zeros(totnodes,1);
Fe = F;
%identify the nodes that are going to be deleted 
del_nodes_2=[1:1:nelemx+1];
del_nodes_4=[(totnodes-nelemy):1:totnodes];

%vector of values of the g_d functions on the nodes that are being deleted
GD_2 = 0;
GD_4 = 0;

for (i=del_nodes_2)
    GD_2 = gd_2( nodes(i,1), nodes(i,2)); 
    Fe = Fe - K(:,i)*GD_2 ;
    
end

for (i=del_nodes_4)
    GD_4 = gd_4( nodes(i,1), nodes(i,2));
    Fe = Fe - (K(:, i))*GD_4 ;
    
end

end

%% 



















