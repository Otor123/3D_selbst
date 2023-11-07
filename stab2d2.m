%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ke,Re,Ae,Ve]=stab2d(e_mat,e_spa,a_var,mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stiffness matrix Ke and righthand side Re for strut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emod = mat(1);   a0 = a_var;   
Wa = mat(3);   
expm = mat(4);   expn = mat(5);   dt   = mat(6); c_rs=mat(7);

indx=[1;3];  ex_mat=e_mat(indx); %x- freiheitsgrad der knoten 
indy=[2;4];  ey_mat=e_mat(indy); %y-freiheitsgrade der knoten
%indz=[]; ez_mat = e_mat(indz); %z-freiheitsgrade der knoten


%Knotenverschiebungen
de=e_spa-e_mat; %e_spa wird verschoben, e_mat - t=0

%Geometrievektor +stablänge // Ursprungsknoten t=0
ax = e_mat(1);
ay = e_mat(2);
bx = e_mat(3);
by = e_mat(4);
dx = ax-bx;
dy = ay-by;
if dx ~= 0 % wenn x verschiebung ist nicht 0
    if dy ~= 0 %Fall dx ~=0 && dy ~= 0!
         phi = atan(dy/dx);
         c = cos(phi);
         s = sin(phi);
    else %Fall dx ~=0 && dy == 0!
         phi = 0;
         c = 1;
         s = 0;
    end
else %Fall dx == 0 && dy ~= 0!
     phi =pi /2; %90 Grad
     c=0;
     s=1;
end
l = sqrt(dx^2+dy^2);
Be= [-c;-s;+c;+s];

%Transformationsmatrix %Werkle 21 S.95

epsi=1/l*Be'*de'; %Dehnung stab %Formel 5.2 Ina Diss S.51
%epsi = 1/l^2 * [dx; dy; -dx; -dy]'*de' 

%update der Stabdicke
var = a_var(1);  %var = Stabdicke
 [var,dres]=updt_2a(epsi,var,mat);

Ae=var;
Ve=Ae*l; %veränderter volumen des stabes

%Normalkräfte an den Knotenkoordinantenkomponenten
Re=(Ae)*emod/l*Be*Be'*de'; %Formel 5.5

% Elementsteifigkeitsmatrix für Stabdicken als Designvariablen
Kee =(Ae)* (emod/l) * Be*Be';
Ke=Kee+(emod/l*Be*Be'*de')*(1/dres)*(c_rs*Ae*emod*1/l^2*Be'*(Be'*de')*dt); %Formel 5.6
%(c_rs*Ae*emod*1/l^2*Be'*(Be'*de')*dt) = c_rs*Ae*emod*1/l^2*1/l * [dx; dy; -dx; -dy]'*(1/l * [dx; dy; -dx; -dy]'*de')*dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%