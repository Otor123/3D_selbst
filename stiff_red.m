function [S]=stiff_red(K,u_pre,ndof)
%Berechnung der reduzierten Steifigkeitsmatrix S_RVE
%alle dofs der Randknoten müssten der ersten Spalte von u_pre entsprechen
rand=u_pre(:,1);
mitte=[1:ndof]';

mitte(rand)=[];

Kpp=K(rand,rand);
Kpr=K(rand,mitte);
Krr=K(mitte,mitte);
Krp=K(mitte,rand);

S=Kpp-Kpr/(Krr)*Krp;






