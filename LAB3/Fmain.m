function [ eL2_global ] = Fmain( Nnodesx, Nnodesy, T1, NT )

%%%%%%%%%%%%%%DATA DEFINITION
q = @(x,y,t) -sin(x-y+t)+(x+y)*cos(x-y+t);
gd_2 = @(x,y,t) +cos(x+t);
gd_3 = @(x,y,t) +cos(1-y+t);
gn_1 = @(x,y,t)-y*sin(-y+t);
gn_4 = @(x,y,t)-x*sin(x-1+t);
in_cond = @(x,y,t) cos(x-y);
ex_sol = @(x,y,t) +cos(x-y+t);

%first and last node of the grid 
x1 = 0;
y1 = 0;
xn = 1;
yn = 1;
%number of nodes (per row, per columns, total)
nnodesx = Nnodesx;
nnodesy = Nnodesy;
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
t1 = T1;
time = linspace(t0,t1,NT);
dt = time(2);

%%%%%%%SOLUTION OF THE PROBLEM 

%I. CREATION OF C MATRIX - time invariant 
C = Cmatrix(nelemx, nelemy, nodes, elements);

%II. INITIAL CONDITION VECTOR
Uini = ICinterp( C, in_cond, totnodes, nodes, elements, t0, nelemx, nelemy);
Uini_proj = ICproj(C, in_cond, totnodes, nodes, elements, t0, nelemx, nelemy);
Uold = Uini;


for t = time
    
    %%%%%%%%%%%%%%STIFFNESS MATRIX K AND LOAD VECTOR F  
    K = zeros (totnodes, totnodes);
    F = zeros (totnodes, 1);
    
    for i=1:totelem
        Ke = LocalstiffK( i, nodes, elements);
        Fe = LocalstiffF(i, nodes, elements, q, t);
        K = K + Ke;
        F = F + Fe;
    end

    %%%%%%%%%%%%%%BOUNDARY CONDITIONS 

    %%%%NEUMANN

    %Gamma vector for boundary 1
    G1 = zeros(totnodes, 1);
    for (i=1:nelemx:totelem)    
        Ge = Boundary_Neum_1(i, nelemx, nelemy, nodes,nnodesx, nnodesy, elements, gn_1, t);  
        G1 = G1 + Ge;    
    end
    %Gamma vector for boundary 4
    G4 = zeros(totnodes, 1);
    for (i=totelem-nelemx+1:1:totelem)    
        Ge = Boundary_Neum_4(i, nelemx, nelemy, nodes,nnodesx, nnodesy, elements, gn_4, t);    
        G4 = G4 + Ge;    
    end
    %Total Gamma vector
    G = G1 + G4;
    %Update the load vector 
    Fnew = F - G;


    %%%%DIRICHLET

    %Add Dirichlet BCs - the dimentions of the system remain the same
    Fdi = Boundary_Dirch_F(Fnew, K, nelemx, nelemy, nodes, gd_2, gd_3, t);
    Kdi = Boundary_Dirch_K(K, nelemx, nelemy, nodes, gd_2, gd_3, t);
    
    %%%%%%%%%%%%%%SOLUTION OF THE SYSTEM
    A = C + dt.*Kdi;
    b = C*Uold + dt.*Fdi;
    
    Unew = pcg(A , b , 1e-12, 1000 );
    Uold = Unew;
    
     

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


%Initial condition in matrix form
% ii = 1;
% Uinimatrix = zeros(nnodesx, nnodesy);
% for i = 1:nnodesy
%     for j=1:nnodesx
%         Uinimatrix(i, j) = Uini(ii);
%         ii = ii+1;
%     end
% end
% ii = 1;
% Uinimatrix_interp = zeros(nnodesx, nnodesy);
% for i = 1:nnodesy
%     for j=1:nnodesx
%         Uinimatrix_interp(i, j) = Uini_interp(ii);
%         ii = ii+1;
%     end
% end

%%%%%%%%%%%%%%REAL SOLUTION 
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

%%%%%%%%%%%%%%PLOTS

figure
subplot(1, 2, 1)
mesh(Ureal)
title('real solution')
subplot(1,2,2)
mesh (Umatrix)
title('numerical solution')
% subplot(1,4,3)
% mesh (Uinimatrix)
% title('initial conditions projected')
% subplot(1,4,4)
% mesh (Uinimatrix_interp)
% title('initial conditions interpolated')



%%%%%%%%%%%%%%ERROR COMPUTATION 
eL2_global = 0;
for (e=1:totelem)
    eL2_local = L2_local_error (e, ex_sol, Uold, nodes, elements, totnodes, time(end));
    eL2_global = eL2_global + eL2_local;
end
eL2_global


end

