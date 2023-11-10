
%Paul Scheuerlein
%Erweiterung der Trab_evo function von 2D in 3D





% %Für 2D:
F = [1.0907, 0; 1.0907, 0 ];

q0_me = [-0.5,-0.5,0.5,0.5; ...
]; %Für jeden IP Koordinaten der Knoten



 randis = [1,2]; %Randknoten
 edof_me = [1 1 2 3 4]; %Beschreibung des einzelknotens
 mat_me = [10000,0.00682274642960738,0.00560000000000000,3,2,0.0100000000000000,133.333333333333]; %Materialpara
 nel_me = 1; %Anzahl der stäbe
 ndof_me = 4; %Anzahl der dof für alle Meso knoten
 a_alt = [0.0068;0.0068;0.0068;0.0068]; 
 e_mat_me=extr_dof(edof_me,q0_me); %jeder knoten in gliabeln KS
 

 %Umrechnung der Deformation in Verschiebungen auf die Randknoten
  [u_pre] = gitter_umrechnung(q0_me(1,:)',F,randis); %Berechnung der Verschiebung eines Elements

%  %3D:
 F = [1.0454, 0, 0; ...
      1.0454, 0, 0; ...
      1.0454, 0, 0];

q0_me = [-0.5,-0.5,-0.5,0.5,0.5,0.5]; 

 randis = [1,2]; %Randknoten
 edof_me = [1 1 2 3 4 5 6]; %Beschreibung des einzelknotens
 mat_me = [10000,0.00682274642960738,0.00560000000000000,3,2,0.0100000000000000,133.333333333333]; %Materialpara
 nel_me = 1; %Anzahl der stäbe
 ndof_me = 6; %Anzahl der dof für alle Meso knoten
 a_alt = [0.0068;0.0068;0.0068;0.0068]; 
 e_mat_me=extr_dof(edof_me,q0_me); %jeder knoten in gliabeln KS

%Umrechnung der Deformation in Verschiebungen auf die Randknoten
  [u_pre] = gitter_umrechnung3d(q0_me(1,:)',F,randis); %Berechnung der Verschiebung eines Elements

  %randis sind die randknoten
  %FE - Berechnung des Mesomodells
  [e_mat_me,a_neu,dof_me,V,P,S_macro,sigma]=...
    trab_evo(q0_me(1,:)',edof_me,u_pre,mat_me,nel_me,ndof_me,a_alt(1,:),F,e_mat_me);