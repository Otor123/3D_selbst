%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ke,Re,Ae,Ve]=stab3d3(e_mat,e_spa,a_var,mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stiffness matrix Ke and righthand side Re for strut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emod = mat(1);   a0 = a_var;   
Wa = mat(3);   
expm = mat(4);   expn = mat(5);   dt   = mat(6); c_rs=mat(7);

indx=[1;4];  ex_mat=e_mat(indx); %x- freiheitsgrad der knoten
indy=[2;5];  ey_mat=e_mat(indy); %y-freiheitsgrade der knoten
indz=[3;6]; ez_mat = e_mat(indz); %z-freiheitsgrade der knoten


%Knotenverschiebungen
de=e_spa-e_mat;

%Geometrievektor +stablänge
ax = e_mat(1);
ay = e_mat(2);
az = e_mat(3);
bx = e_mat(4);
by = e_mat(5);
bz = e_mat(6);

dx = bx-ax;
dy = by-ay;
dz = bz-az;


% if dx~= 0 && dy~= 0 && dz~=0
% 
%     alpha = atan(dy/dx); 
%     %mit alpha weil dz bei alpha und beta
%     %beta = atan(dy/dx); 
%     c = cos(alpha);
%     s = sin(alpha);
% 
% elseif (dx~=0 && dy~=0 && dz==0)
% 
% 
%     alpha = atan(dy/dx); 
%     c = cos(alpha);
%     s = sin(alpha);
% 
% 
% elseif (dx~=0 && dy==0 && dz~=0)
% 
%     beta = atan(dx/dz); 
%     c = cos(beta);
%     s = sin(beta);
% 
% elseif (dx~=0 && dy==0 && dz==0)
% 
%     beta = pi/2; 
%     c = cos(beta);
%     s = sin(beta);
% 
% elseif (dx==0 && dy~=0 && dz~=0)    
% 
%     alpha = atan(dy/dx); 
%     c = cos(alpha);
%     s = sin(alpha);
% 
% elseif (dx==0 && dy~=0 && dz==0)  %NOCH NICHT
% 
%     alpha = atan(dy/dx); 
%     c = cos(alpha);
%     s = sin(alpha);    
% 
% elseif (dx==0 && dy==0 && dz~=0)      %NOCH NICHT
% 
%     alpha = atan(dy/dx); 
%     c = cos(alpha);
%     s = sin(alpha);    
% end


l = sqrt(dx^2+dy^2+dz^2);

de_lok = 1/l * [dx dy dz 0 0 0; 0 0 0 dx dy dz] * de';

dl = de_lok(2) - de_lok(1);
%Be= [-c;-s;+c;+s];

epsi=1/l*dl; %Be'*de'; %Dehnung stab %Formel 5.2
 %epsi = 1/l^2 * [-dx; -dy; -dz; dx; dy; dz]'*de' 
%update der Stabdicke
var = a_var(1);  
[var,dres]=updt_2a(epsi,var,mat);

Ae=var;
Ve=Ae*l;

%Elementsteifigkeitsmatrix
Kee = (emod*Ae/l^3)* [dx^2 dx*dy dx*dz -dx^2 -dx*dy -dx*dz; ...
                     dx*dy dy^2 dy*dz -dx*dy -dy^2 -dy*dz; ...
                     dx*dz dy*dz dz^2 -dx*dz -dy*dz -dz^2; ...
                     -dx^2 -dx*dy -dx*dz dx^2 dx*dy dx*dz; ...
                     -dx*dy -dy^2 -dy*dz dx*dy dy^2 dy*dz; ...
                     -dx*dz -dy*dz -dz^2 dx*dz dy*dz dz^2];


%Normalkräfte an den Knotenkoordinantenkomponenten
Re = Kee * de';
%Re=((Ae)*emod/l)*Be*Be'*de';
%Re=((Ae)*emod/l^3)*[-dx -dy -dz dx dy dz]*[-dx -dy -dz dx dy dz]'*de'
% Elementsteifigkeitsmatrix für Stabdicken als Designvariablen
%Kee =(Ae)* (emod/l) * Be*Be';
%Ke=Kee+(emod/l*Be*Be'*de')*(1/dres)*(c_rs*Ae*emod*1/l^2*Be'*(Be'*de')*dt);

Ke = Kee+(Re/Ae)*(1/dres)*(c_rs*Ae*emod*1/l^2* 1/l * [-dx; -dy; -dz; dx; dy; dz]'*(dl)*dt);

%(Be'*de') == dl vgl Z 106 u. 86
%Be == Kee/(Ae)*emod/l vgl Z 106 mit Kee == ((Ae)*emod/l)*Be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%