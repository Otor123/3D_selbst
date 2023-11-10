%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q0,edof,emat,nel,node,ndof,V_qa,knoten,staebe,a0,mat,randis] = rve_ini(art,mat,zuf,anz,iny)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Erstellung des Initialen Stabmodells und der Meso-Parameter
%art: Gitterart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% material parameters for strut evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rho0=mat(3);  %homogenisierte Anfangsdichte entspricht macro Anfangsdichte

emod=mat(1);  %Entspricht dem solid_emodul wie in Macro
Wa=mat(4)*4;%fem*7  %Laut jacob etwa mal faktor 7 o 6.7  femur:*3 square=mat(4)   neu: square*4
expm=mat(5);  %wie in macro
expn=mat(6);  %wie in macro
dt=mat(7);    %wie in macro
c_rs=mat(8)/1.5; %fem/7.5    %laut jacobb etwa durch faktor 20 kleiner   femur:1/20 neu:  square=mat(8)/4 neu/1.5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate node points q0 and connectivity edof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Größe des repräsentativen Volumenelements
xbox(1) = 0;
xbox(2) = 1;           %ca 4 mal kleiner als macro element
ybox(1) = 0;
ybox(2) = 1;
zbox(1) = 0;
zbox(2) = 1;

%anzahln der Knoten in eine Richtung
if nargin >3
   nx=anz;
   ny=anz;
elseif nargin > 4
    nx=anz;
    ny=iny;
else
    nx =  10;
    ny =  10;
end

%Fläche des RVE
V_qa=(xbox(2)-xbox(1))*(ybox(2)-ybox(1)*(zbox(2)-zbox(1)));


%Gittererstellung (nach Grossberger)
[knoten,staebe]=gridgeneration(xbox(1),xbox(2),ybox(1),ybox(2),zbox(1),zbox(2), nx,ny,nz,art,zuf);
mittexy=(xbox(2)-xbox(1))/2;
knoten=knoten-mittexy;
xbox=xbox-mittexy;
ybox=ybox-mittexy;

%Umrechnung der Knotenkoord. (nodes,2) in einen Spaltenvektor (2*nodes,1)
for i=1:size(knoten,1)
    q0(2*i-1,1)=knoten(i,1);
    q0(2*i,1)=knoten(i,2);
end

%q0 - enthaält alle x,y werte nacheinander

%Zuordnung der Freiheitsgrade zu Elemente
for i=1:size(staebe,1)
    edof(i,1)=i;
    edof(i,2)=staebe(i,1)*2-1; %anfangsknote - nr freiheitsgrad x-richtung
    edof(i,3)=staebe(i,1)*2; %anfangsnkoten - nr freiheitsgrad y-richtung
    edof(i,4)=staebe(i,2)*2-1; %endknoten - nr freiheitsgrad x-richtung
    edof(i,5)=staebe(i,2)*2; %endknoten - nr freiheitsgrad y-richtung
end

%Wichtige Größen für den Algo.
[nel, sizen] = size(edof);
[ndof,sizen] = size(q0);
 for ie=1:nel emat(ie) = 1; end; %Zeilenvektor mit 1en
 node        = ndof/2;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initiale Trabekeldicke bestimmen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%stablängen Berechnung
[sn,sr] = size(staebe);
laengen = zeros(sn,1) ;
l_ges=0;
for s = 1:sn
    a = staebe(s,1); %anfangsknoten
    b = staebe(s,2); %endknoten
    ax = knoten(a,1); %x-pos. anfangsknoten i
    ay = knoten(a,2); %y-pos. anfangsknoten i
    bx = knoten(b,1); %x-pos. endknoten i
    by = knoten(b,2); %y-pos. endknoten i
    dx = ax-bx; 
    dy = ay-by;
    laengen(s) = sqrt(dx^2+dy^2);
    l_ges=l_ges+laengen(s);
end

rhos=2;   %solid dichte der Stäbe

a0=rho0/rhos*V_qa/l_ges;            %Berechnung des initialen Stabquerschnitts

mat = [emod,a0,Wa,expm,expn,dt,c_rs]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Randknotenbestimmung

z=0;
for k=1:size(knoten,1)                      %randknoten herausfinden
    if knoten(k,1)==xbox(1)
        z=z+1;
        randknoten(z,:)=knoten(k,:);        %knoten befindet sich am linken rand
        randis(z)=k;
    elseif knoten(k,1)==xbox(2)
        z=z+1;
        randknoten(z,:)=knoten(k,:);        %knoten befindet sich am rechten rand
        randis(z)=k;
    elseif knoten(k,2)==ybox(1)
        z=z+1;
        randknoten(z,:)=knoten(k,:);        %knoten befindet sich am unteren rand
        randis(z)=k;
    elseif knoten(k,2)==ybox(2)
        z=z+1;
        randknoten(z,:)=knoten(k,:);        %knoten befindet sich am oberen rand
        randis(z)=k;
    else
    end
end