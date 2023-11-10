

%Erweiterung rve_ini für 3D FEM
%Paul Scheuerlein
%09.11.2023


rve_ini(2,mat,0,9);
mat = [10000,0.330000000000000,0.400000000000000,0.001400000000000,3,2,0.010000000000000,200];

[q0,edof,emat,nel,node,ndof,V_sq,knoten,staebe,a0,mat,randis] = rve_ini(2,mat,0,9);