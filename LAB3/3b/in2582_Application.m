%This code is meant to solve parabolic problems with boundary codnitions =0, and a point-source function in the center of the domain
%It combines a finite element method with
%backard Euler for the time computations. The code is able to handle sqare
%and rectangular elements, since the number of nodes on x and y axis can be
%different.
% { line 30 - Diffusion coeffcient
% { line 90 - Forcing function
% { line 84 - Dirichlet boundary values
% { line 44 - Reference Finite element
% { line 26 - Initial conditions
% { line 60 - Projection of initial data
% { line 88, 89 - Point source



clc 
clear all 


%%%%%%%%%%%%%%DATA DEFINITION
q = @(x,y,t) 0; 
gd_2 = @(x,y,t) 0;
gd_3 = @(x,y,t) 0;
gd_1 = @(x,y,t) 0;
gd_4 = @(x,y,t) 0;
in_cond = @(x,y) 0;
ex_sol = @(x,y,t) +cos(x-y+t);
k = @(x,y) 1*(x<=0.5) + 10*(x>0.5);
%first and last node of the grid 
x1 = 0;
y1 = 0;
xn = 1;
yn = 1;
%number of nodes (per row, per columns, total)
nnodesx = 11;
nnodesy = 11;
totnodes = nnodesx*nnodesy;
%number of elements (per row, per column, total)
nelemx = nnodesx - 1;
nelemy = nnodesy - 1;
totelem = nelemx*nelemy;

%vector of nodes with their corrispective coordinates
[ nodes ] = RNODE( nnodesx, nnodesy, x1, y1, xn, yn );
%matrix in which elements are linked to their corrispective nodes 
[ elements ] = RELEM( totelem, nelemx );

%time discretization
t0 = 0;
t1 = 1;
time = linspace(t0,t1,11);
dt = time(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLUTION OF THE PROBLEM 

%%%%%%%%%%%%%% I. CREATION OF C MATRIX - time invariant 
C = Cmatrix(nelemx, nelemy, nodes, elements);

%%%%%%%%%%%%%% II. INITIAL CONDITION VECTOR
Uini = ICinterp( C, in_cond, totnodes, nodes, elements);
Uold = Uini;

kk = floor(nnodesx/2+1);
U_middle = zeros(nnodesy, 4);
U_middle(:,1) = Uold( kk:nnodesx:(totnodes-kk+1));

for t = time
    
    %%%%%%%%%%%%%% III. STIFFNESS MATRIX K AND LOAD VECTOR F  
    K = zeros (totnodes, totnodes);
    F = zeros (totnodes, 1);
    
    for i=1:totelem
        Ke = LocalstiffK( i, nodes, elements, k);
        Fe = LocalstiffF(i, nodes, elements, q, t);
        K = K + Ke;
        F = F + Fe;
    end

    %%%%%%%%%%%%%% IV. BOUNDARY CONDITIONS 
    
    %%%%DIRICHLET

    %Add Dirichlet BCs - the dimentions of the system remain the same -
    %here the source founction has been inserted directly into the
    %functions - It is handle like a dirichlet condition in the central
    %node of the mesh. 
    Fdi = Boundary_Dirch_F(F, K, nelemx, nelemy, nodes, gd_1, gd_2, gd_3, gd_4, t);
    Kdi = Boundary_Dirch_K(K, nelemx, nelemy);
    
    %%%%%%%%%%%%%% V.SOLUTION OF THE SYSTEM
    A = C + dt.*Kdi;
    b = C*Uold + dt.*Fdi;
    
    Unew = pcg(A , b , 1e-12, 1000 );
    Uold = Unew;
    
    %Memorization of value on x_2 = 0.5
    if (t==0.005)
        U_middle(:,2) = Unew( kk:nnodesx:(totnodes-kk+1));
    elseif(t==0.01)
        U_middle(:,3) = Unew( kk:nnodesx:(totnodes-kk+1));
    elseif(t==0.02)
        U_middle(:,4) = Unew( kk:nnodesx:(totnodes-kk+1));
    end
    
     

end


%Solution in matrix form
ii = 1;
Umatrix = zeros(nnodesx, nnodesy);
for i = 1:nnodesy
    for j=1:nnodesx
        Umatrix(i, j) = Uold(ii);
        ii = ii+1;
    end
end

% %Initial condition in matrix form
% ii = 1;
% Uinimatrix = zeros(nnodesx, nnodesy);
% for i = 1:nnodesy
%     for j=1:nnodesx
%         Uinimatrix(i, j) = Uini(ii);
%         ii = ii+1;
%     end
% end

%%%%%%%%%%%%%% VI. REAL SOLUTION 
%matrix of the real solution values
Ureal = zeros(nnodesx, nnodesy);
ureal = zeros(totnodes, 1);
ii=1;
for (i = 1:nnodesy)
    for(j=1:nnodesx)        
        Ureal(i,j) = ex_sol(nodes(ii, 1), nodes(ii,2), time(end));
        ii = ii+1;
    end
end
% for (i = 1:totnodes)
%     ureal(i) = ex_sol(nodes(i, 1), nodes(i,2), time(end));
% end

%%%%%%%%%%%%%% VII. PLOTS

% figure
% subplot(1, 2, 1)
% mesh(Ureal)
% title('real solution')
% subplot(1,2,2)
imagesc(Umatrix)
title('numerical solution')
xlabel('x')
ylabel('y')
% subplot(1,3,3)
% mesh (Uinimatrix)
% title('initial conditions')



%%%%%%%%%%%%%% VII. ERROR COMPUTATION 
eL2_global = 0;
for (e=1:totelem)
    eL2_local = L2_local_error (e, ex_sol, Uold, nodes, elements, totnodes, t);
    eL2_global = eL2_global + eL2_local;
end
eL2_global

%%
% x = linspace(0,1, 21);
% plot(x, U_middle(:,1))
% hold on
% plot(x, U_middle(:,2))
% hold on
% plot(x, U_middle(:,3))
% hold on
% plot(x, U_middle(:,4))
% legend ('t=0', 't=0.005', '0.01', '0.02')
%%
function [ nodes ] = RNODE( nnodesx, nnodesy, x1, y1, xn, yn )

nelemx = nnodesx - 1;
nelemy = nnodesy - 1;
%definition of increments = side length of each element 
dx = (xn-x1)/nelemx;
dy = (yn-y1)/nelemy;
totnodes = nnodesx*nnodesy;
%vector of nodes coordinates
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
function [ C ] = Cmatrix( nelemx, nelemy, nodes, elements)
%this fucntion creates the C matrix - time invariant 

totnodes = (nelemx+1)*(nelemy+1);
totelem = nelemx*nelemy;

C = zeros(totnodes, totnodes);

for e = 1:totelem
    
        Ce = zeros(totnodes, totnodes);
        cur_nodes = nodes(elements(e, :)',:);
        
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

        %weight functions
        w1 = 1;
        w2 = 1;
        w3 = 1;
        w4 = 1;
        
        %basis functions at integration points 
        PHI = zeros(4,4);
        DXDE = zeros(4,1);
        DXDN = zeros(4,1);
        DYDE = zeros(4,1);
        DYDN = zeros(4,1);
        J = zeros(4,1);
        
        for i=1:4
            for j=1:4
            PHI(i,j) = master_shape(j, int_points(i,1), int_points(i,2));
            DXDE(i) = DXDE(i) + cur_nodes(j,1)*master_shape_deriv_x(j,int_points(i,1),int_points(i,2));
            DXDN(i) = DXDN(i) + cur_nodes(j,1)*master_shape_deriv_y(j,int_points(i,1),int_points(i,2)); 
            DYDE(i) = DYDE(i) + cur_nodes(j,2)*master_shape_deriv_x(j,int_points(i,1),int_points(i,2));
            DYDN(i) = DYDN(i) + cur_nodes(j,2)*master_shape_deriv_y(j,int_points(i,1),int_points(i,2)); 
            end
        end
        
        for i=1:4
            J(i) = DXDE(i)*DYDN(i)-DXDN(i)*DYDE(i);
        end
        
        w1 = w1*J(1);
        w2 = w2*J(2);
        w3 = w3*J(3);
        w4 = w4*J(4);
        
        ii = 1;             
        for i = elements(e, :)
            jj = 1;
            for j = elements(e, :)
                
                %this is TO CHECK!! 
                k1 = PHI(1,ii)*PHI(1,jj);
                k2 = PHI(2,ii)*PHI(2,jj);
                k3 = PHI(3,ii)*PHI(3,jj);
                k4 = PHI(4,ii)*PHI(4,jj);

                Ce(i,j) = w1*k1+w2*k2+w3*k3+w4*k4 ;
                
                jj = jj+1;
            end
            ii = ii+1;
        end
        
     C = C + Ce;   
        
end

end

%%
function [ Uold ] = ICinterp( C, in_cond, totnodes, nodes, elements)

Uold = zeros(totnodes, 1);

for (i = 1:totnodes)
    Uold(i) = in_cond( nodes(i,1), nodes(i,2) );
end

end

%%
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
%%
function [ Fe ] = LocalstiffF( e, nodes, elements, q, t )

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
    F(i) = q(int_points(i,1), int_points(i,2), t);
    K(i*2-1:i*2, :) = [Y(i), 0 ; 0, X(i)];
end

%compute the function in the 4 integration points and then sum them. do the
%same fot the load vector. 

w1 = w1*J(1);
w2 = w2*J(2);
w3 = w3*J(3);
w4 = w4*J(4);

ii=1;

for (i=elements(e, :))
    
      %i is the node
      %ii is the index of the basis function (1:4)
      f1 =  q(X(1), Y(2), t)*PHI(1, ii);
      f2 =  q(X(2), Y(2), t)*PHI(2, ii);
      f3 =  q(X(3), Y(3), t)*PHI(3, ii);
      f4 =  q(X(4), Y(4), t)*PHI(4, ii);      
      Fe(i) = w1*f1+w2*f2+w3*f3+w4*f4;
      ii = ii+1;
   
end




end
%%

function [ K ] = Boundary_Dirch_K(K, nelemx, nelemy)

%identify the nodes that are going to be deleted 
nnodesx = nelemx+1;
nnodesy = nelemy+1;
totnodes = (nelemx+1)*(nelemy+1);
del_nodes_2 = [1:1:nelemx+1];
del_nodes_3 = [2*(nelemx+1) : nelemx+1 : totnodes];
del_nodes_1 = [nnodesx+1:nnodesx:totnodes-nnodesx+1];
del_nodes_4 = [totnodes-nnodesx+1:1:totnodes];

for (i=del_nodes_2)
    K(i,:) = zeros(1, totnodes);
    K(:,i) = zeros(totnodes,1);
    K(i,i) = 1;
end

for (i=del_nodes_3)
    K(i,:) = zeros(1, totnodes);
    K(:,i) = zeros(totnodes,1);
    K(i,i) = 1;
end

for (i=del_nodes_1)
    K(i,:) = zeros(1, totnodes);
    K(:,i) = zeros(totnodes,1);
    K(i,i) = 1;
end

for (i=del_nodes_4)
    K(i,:) = zeros(1, totnodes);
    K(:,i) = zeros(totnodes,1);
    K(i,i) = 1;
end

%centered source
k = floor(totnodes/2)+1;
K(k,:) = zeros(1, totnodes);
K(:,k) = zeros(totnodes,1);
K(k,k) = 1;



end
%%

function [ Fe ] = Boundary_Dirch_F(F, K, nelemx, nelemy, nodes, gd_1, gd_2, gd_3, gd_4, t)
%this function modifies the load vector adding the Dirichlet boudnary
%conditions

nnodesx=nelemx+1;
totnodes = (nelemx+1)*(nelemy+1);
totelem = nelemx+nelemy;
Fe = zeros(totnodes,1);
Fe = F;

%identify the nodes that are going to be deleted 
del_nodes_2 = [1:1:nelemx+1];
del_nodes_3 = [2*(nelemx+1) : nelemx+1 : totnodes];
del_nodes_1 = [nnodesx+1:nnodesx:totnodes-nnodesx+1];
del_nodes_4 = [totnodes-nnodesx+1:1:totnodes];

%vector of values of the g_d functions on the nodes that are being deleted

for (i=del_nodes_1)
    GD_1 = gd_1( nodes(i,1), nodes(i,2), t); 
    Fe = Fe - K(:,i)*GD_1 ;
    
end

for (i=del_nodes_2)
    GD_2 = gd_2( nodes(i,1), nodes(i,2), t); 
    Fe = Fe - K(:,i)*GD_2 ;
    
end

for (i=del_nodes_3)
    GD_3 = gd_3( nodes(i,1), nodes(i,2), t);
    Fe = Fe - (K(:, i))*GD_3 ;
    
end
for (i=del_nodes_4)
    GD_4 = gd_4( nodes(i,1), nodes(i,2), t); 
    Fe = Fe - K(:,i)*GD_4 ;
    
end

% Now I can substitue the right values 
for (i = del_nodes_1)
    Fe(i) = gd_1( nodes(i,1), nodes(i,2), t);
end
for (i = del_nodes_2)
    Fe(i) = gd_2( nodes(i,1), nodes(i,2), t);
end
for (i = del_nodes_3)
    Fe(i) = gd_3( nodes(i,1), nodes(i,2), t);
end
for (i = del_nodes_4)
    Fe(i) = gd_4( nodes(i,1), nodes(i,2), t);
end

%centered source 
k = floor(totnodes/2)+1;
if (t<=0.01)
    GD = 1+100*t; 
    Fe = Fe - K(:,k)*GD ;
    Fe(k) = 1+100*t;
elseif(t>0.01)
    GD = 2; 
    Fe = Fe - K(:,k)*GD ;
    Fe(k) = 2;
end


end
%%
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




