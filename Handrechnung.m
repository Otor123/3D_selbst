clc
clear all

%Paul Scheuerlein
%Erweiterung des Stab2d2 zu Stab3d3
%Quelle: Werkle, 21
ende = 70;
N = +1000; %+ Zug, - Druck
start_F1 = -N; %immer Zeile 22 beachten
start_F2 = 1000;
dF = N/5;
dF1 =  -dF;
dF2 =  dF;
F1 = start_F1;
F2 = start_F2;


A_plot  = zeros(1,ende);
F_plot1 = zeros(1,ende);
F_plot2 = zeros(1,ende);


for m = 1:ende
%F1 und F2 stellt die lokalen Stabkröfte dar
F1 =  F1 + (m-1)*dF1;
F2 =  F2 + (m-1)*dF2;
%3D - Berechnung der Globalen kräfte: 
F1i = sqrt(1/3*(F1^2));
F2i = sqrt(1/3*(F2^2));
%2D - Berechnung der Globalen kräfte:
F1i2D = sqrt(1/2*(F1^2));
F2i2D = sqrt(1/2*(F2^2));
%E-Modul
E = 10000;
%Querschnitt
A = 0.0068;
%Koordinaten von Knoten 1 und Knoten 2
k1 = [-0.5,-0.5,-0.5];
k2 = [0.5,0.5,0.5];
%Berechnung der delta längen
lx = k2(1)-k1(1);
ly = k2(2)-k1(2);
lz = k2(3)-k1(3);
%Berechnung der Stablängen für 3 und 2D
l = sqrt(lx*lx+ly*ly+lz*lz);
l_2d = sqrt(lx*lx+ly*ly);

%Berechnung Elementsteifigkeitsmatrix des Fachwerkstabes 1D
% Ke1D = E*A/l * [1 -1];
% u1D = F1\Ke1D;

Ke1D_alter = E*A/l * [1 -1; -1 1];
F1D_alter = [F1; F2];
u1D_alter = F1D_alter \ Ke1D_alter;



%Berechnung Elementsteifigkeitsmatrix des Fachwerkstabes 2D

%stimmt an sich // zu viel aufwand zum vergleich


% KE2D_adv = (E*A/(l_2d^3)) * [lx^2 lx*ly -lx^2 -lx*ly; ...
%                             lx*ly ly^2 -lx*ly -ly^2; ...
%                             -lx^2 -lx*ly lx^2 lx*ly; ...
%                             -lx*ly -ly^2 lx*ly ly^2];
% 
%  
% 
% F2D_adv = [-F1i2D; -F1i2D; F2i2D; F2i2D];
% 
% u2D_adv = F2D_adv \ KE2D_adv;


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
    %global:
U3D = Fe \ Ke;
    %lokal:
u_lok_3D = 1/l * [lx, ly, lz, 0, 0, 0; 0, 0, 0, lx, ly, lz] * U3D';

% u_3D_lok(1) = -sqrt((U3D(1)^2)*3);
% u_3D_lok(2) = sqrt(U3D(4)^2*3);

%Kontrolle der Kräfte und Verschiebungen im CommandWindow
fprintf('\nVerschiebung: 3D-lokal     3D-global (x,y,z)             1D-lokal  \n' )
fprintf('Node1: %15.5f  %12.5f %5.5f %5.5f  %9.5f      \n' , u_lok_3D(1), U3D(1), U3D(2), U3D(3), u1D_alter(1) )
fprintf('Node2: %15.5f  %12.5f %8.5f %8.5f  %9.5f      \n' , u_lok_3D(2), U3D(4), U3D(5), U3D(6), u1D_alter(2) )
fprintf('\nKraft lokal: F1_lok = %5.2f F2_lok = %5.2f \n', F1, F2 );


%Für Plot
verschiebung(m,1) = u_lok_3D(1);
verschiebung(m,2) = u_lok_3D(2);
%e_mat_me=extr_dof(edof_me,q0_me);  %nlin_mm Z.9

%Input für stab3d3
edof_3d = [1 1 2 3 4 5 6 ];
q0_3d = [k1, k2];
e_spa_3d = q0_3d + U3D;
e_mat_3d=extr_dof(edof_3d,q0_3d); 
a_var_3d(1) = 0.0068;
e_par_3d = [10000, 0.00682274642960738,0.00560000000000000,3,2,0.0100000000000000,133.333333333333];
%stab3d3 aufruf 
[Ke,Re,Ae] = stab3d3(e_mat_3d(1,:),e_spa_3d(1,:),a_var_3d(1,:),e_par_3d(1,:));

%Input für stab2d2
edof_2d = [1 1 2 3 4];
q0_2d = [k1(1:2), k2(1:2)];
e_spa_3d(3) = [];
e_spa_3d(5) = [];
e_spa_2d = e_spa_3d;
e_mat_2d=extr_dof(edof_2d,q0_2d); 
a_var_2d(1) = 0.0068;
e_par_2d = [10000, 0.00682274642960738,0.00560000000000000,3,2,0.0100000000000000,133.333333333333];
%stab2d2 aufruf 
[Ke_2d,Re_2d,Ae_2d] = stab2d2(e_mat_2d(1,:),e_spa_2d(1,:),a_var_2d(1,:),e_par_2d(1,:));


fprintf('\nneuer Querschnitt 2D [m^2]: %5.5f \n' , Ae_2d )
fprintf('neuer Querschnitt 3D [m^2]: %5.5f \n' , Ae )
A_plot(m) = Ae;
A_plot_2D(m) = Ae_2d;
F_plot1(m) = F1;
F_plot2(m) = F2;
zyklen(m) = m;


end

%Plotten der Kräfte, Querschnitt und Verschiebung

figure(1) %Kräfte
plot(zyklen, F_plot1, zyklen, F_plot2);
title('Kräfte F1 und F2 lokal am Stab')
legend('F1-lok','F2-lok')
xlabel('Last Zyklen')
ylabel('Kraft [N]')

figure(2) %Querschnitt
plot(zyklen, A_plot, zyklen, A_plot_2D);
title('Querschnittsveränderung des Stabes 2D und 3D')
xlabel('Last Zyklen')
ylabel('Querschnitt [m²]')

figure(3) %Verschiebung
plot(zyklen, verschiebung(:,1), zyklen, verschiebung(:,2));
title('Verschiebung des Stabes nur 3D')
xlabel('Last Zyklen')
ylabel('Verschiebung [m]')