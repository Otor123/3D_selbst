clc
clear all

%Paul Scheuerlein
%Quelle: Werkle, 21




F1 =  500;
F2 =  500;
%3D:
F1i = sqrt(1/3*(F1^2));
F2i = sqrt(1/3*(F2^2));
%2D
F1i2D = sqrt(1/2*(F1^2));
F2i2D = sqrt(1/2*(F2^2));

E = 10000;

A = 0.0068;
k1 = [-0.5,-0.5,-0.5];
k2 = [0.5,0.5,0.5];

lx = k2(1)-k1(1);
ly = k2(2)-k1(2);
lz = k2(3)-k1(3);

l = sqrt(lx*lx+ly*ly+lz*lz);
l_2d = sqrt(lx*lx+ly*ly);

%Berechnung Elementsteifigkeitsmatrix des Fachwerkstabes 1D

Ke1D = E*A/l * [-1 1];
u1D = F1\Ke1D;

Ke1D_alter = E*A/l * [1 -1; -1 1];
F1D_alter = [-F1; F2];
u1D_alter = F1D_alter \ Ke1D_alter;



%Berechnung Elementsteifigkeitsmatrix des Fachwerkstabes 2D

%stimmt an sich // zu viel aufwand zum vergleich
% 
% 
KE2D_adv = E*A/l_2d^3 * [lx^2 lx*ly -lx^2 -lx*ly; ...
                      lx*ly ly^2 -lx*ly -ly^2; ...
                      -lx^2 -lx*ly lx^2 lx*ly; ...
                      -lx*ly -ly^2 lx*ly ly^2];

 

F2D_adv = [-F1i2D; 0; F2i2D; 0];

u2D_adv = F2D_adv \ KE2D_adv;


%%Berechnung Elementsteifigkeitsmatrix des Fachwerkstabes 3D

Ke = E*A/(l^3)*...
    [lx*lx,lx*ly,lx*lz,-lx^2,-lx*ly,-lx*lz;...
    lx*ly,ly^2,ly*lz,-lx*ly,-ly^2,-ly*lz;...
    lx*lz,ly*lz,lz^2,-lx*lz,-ly*lz,-lz^2;...
    -lx^2,-lx*ly,-lx*lz,lx^2,lx*ly,lx*lz;...
    -lx*ly,-ly^2,-ly*lz,lx*ly,ly^2,ly*lz;...
    -lx*lz,-ly*lz,-lz^2,lx*lz,ly*lz,lz^2];

Fe = zeros(6,1);



for i = 1:length(Fe)

    if i <= length(Fe)/2
    Fe(i) = -F1i;

    else

    Fe(i) = F2i;
    end
end

%Verschiebung;
U3D = Fe \ Ke;

disp(U3D);

u1_3D_glo = -sqrt((U3D(1)^2)*3);
u2_3D_glo = sqrt(U3D(4)^2*3);



%%FEM NOTWEWNDIG

%e_mat_me=extr_dof(edof_me,q0_me);  %nlin_mm Z.9


edof_ha = [1 1 2 3 4 5 6 ];
q0_ha = [-0.5, -0.5, -0.5, 0.5 ,0.5 0.5];

e_spa_ha = q0_ha + U3D;

e_mat_ha=extr_dof(edof_ha,q0_ha); 

a_var_ha(1) = 0.0086;

e_par_ha = [10000, 0.00682274642960738,0.00560000000000000,3,2,0.0100000000000000,133.333333333333];

%e_spa_ha = 

%nun muss stab2d aufruf folgen

[Ke,Re,Ae] = stab3d3(e_mat_ha(1,:),e_spa_ha(1,:),a_var_ha(1,:),e_par_ha(1,:));
 


