function [knoten,staebe]= gridgeneration3d(gx0,gx1,gy0,gy1,gz0,gz1,nelx,nely,nelz,art,zuf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function - generation of nodes and trabeculae
%
% Date: 2019-10-02
% Worker: AP
% Modified: AG 2021-03-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 10
    art = 0;
end

% determinition of the nodes (coordinates)
% varx= (gx1-gx0)/30.;  % factor for variation of coordinates in x direction
% vary= (gy1-gy0)/10.;  % factor for variation of coordinates in y direction
if zuf>0
    varx= (gx1-gx0)/nelx/4.;  % factor for variation of coordinates in x direction
    vary= (gy1-gy0)/nely/3.;  % factor for variation of coordinates in y direction
    varz= (gz1-gz0)/nelz/3.;
else
    varx=0.;
    vary=0.;
    varz=0.;
end


for n = 1:nelz
jj=0;
j=1;
z=0;
     for i=1:(2*(nelx-2)+2*(nely))
        if i <= (nely+1)                                  % left boundary
            bound(n,i)=i;              
        end
        if i > (nely+1) && i < (2*(nelx-2)+2*(nely))-nely % lower and upper boundary
            bound(n,i)=i + j*(nely-2);               
            z=z+1;
            if z==2
                z=0;
                j=j+1;
            end
        end
        if i >= (2*(nelx-2)+2*(nely))-nely                  % right boundary
            bound(n,i)=(nelx)*(nely)-nely+jj;           
            jj=jj+1;
        end
    end
end

count = 0;
for n = 1:nelz
    bound(n, :) = bound (n, :) + (nelx*nely)*count;
    count = count +1;
end



% Matrix in Vektor umwandeln
bound = reshape(bound.', 1, []);


% for i=1:(2*(nelx-2)+2*(nely))
% bound(i)
% end
fini=1;

%Hierbei erfolgt x,y koordinate für jeden knoten
knoten = zeros (nelx*nely*nelz,3) ;
for elz = 1: nelz
    for elx = 1: nelx
        for ely = 1: nely
            k = ely+nely*(elx -1)+(nelx*nely)*(elz -1); % node number
            b=0;
            for i=fini:(2*(nelx-2)+2*(nely)) %Es erfolgt abfrage ob Randknoten
                if k==bound(i)
                    b=1;   % boundary node
                    fini=i;
                end
            end
            if b==1 % boundary node
                % equidistributed nodes per direction, esp for boundary
                knoten(k,1) = (elx-1)*(gx1-gx0)/(nelx-1); % x-coordinate
                knoten(k,2) = (ely-1)*(gy1-gy0)/(nely-1); % y-coordinate
                knoten(k,3) = (elz-1)*(gz1-gz0)/(nelz-1); % z-coordinate
            else % inner node
                % node coordinates with a random excentricity,esp for inner nodes
                if ismember(art,[0,1,2,4,5,6]) %sollte was mit zufälliger poistion zu tun haben
                    knoten(k,1) = ( (elx-1)*(gx1-gx0)/(nelx-1)) + (-1 + 2.*rand(1))*varx ; % x-coordinate
                    knoten(k,2) = ( (ely-1)*(gy1-gy0)/(nely-1)) + (-1 + 2.*rand(1))*vary ; % y-coordinate
                    knoten(k,3) = ( (elz-1)*(gz1-gz0)/(nelz-1)) + (-1 + 2.*rand(1))*varz ; % z-coordinate
                else %random excentricity for hexagon structures are applied later
                    knoten(k,1) = (elx-1)*(gx1-gx0)/(nelx-1); % x-coordinate
                    knoten(k,2) = (ely-1)*(gy1-gy0)/(nely-1); % y-coordinate
                    knoten(k,3) = (elz-1)*(gz1-gz0)/(nelz-1); % y-coordinate
                end
            end %end of if
        end %end of ely
    end %end of elx
end %end of elz

% Plotten der Knoten
figure;
plot3(knoten(:,1), knoten(:,2), knoten(:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
grid on;
xlabel('X-Achse');
ylabel('Y-Achse');
zlabel('Z-Achse');
title('Knoten');

%Zusätzliche Mitten-Knoten für 6 Ecke
if art == 9 || art == 10
   for elx = 1 : fix((nelx-1)/2 )
      for ely = 1 : nely
          if mod(ely,2)==0 && mod(elx,2)~=0
              k = k + 1;
              knoten(k,1) = 1.5* ((gx1-gx0)/(nelx-1)) + (elx-1) * 2 * ((gx1-gx0)/(nelx-1)) ; % x-coordinate
              knoten(k,2) =  ( (ely-1)*(gy1-gy0)/(nely-1)) ; % y-coordinate
          
          elseif mod(ely,2)~=0 && mod(elx,2)==0
              k = k + 1;
              knoten(k,1) = 1.5* ((gx1-gx0)/(nelx-1))  + (elx-1) * 2 * ((gx1-gx0)/(nelx-1)); % x-coordinate
              knoten(k,2) =  ( (ely-1)*(gy1-gy0)/(nely-1)) ; % y-coordinate
          end
      end %ely
   end %elx
end






% determinition of the trabeculae
staebe = zeros((nelx*nely)^2,2);
sn = 1;
%fprintf('type of grid generation: art = %d \n',art);
if art == 2
    % art == 2 : connections only between neighboured nodes without
    % crossings
for elz = 1:nelz
    for elx = 1:nelx
        for ely = 1:nely
            k1 = ely+nely*(elx-1)+(nelx*nely)*(elz -1); % node number
            

            if ( mod(elx,2) == mod(1+nely-ely,2) ) %wozu??
                ii = [0:1];
                jj = [-1:1];
                k = 0;
            else
                ii = [0:1];
                jj = [0:1];
                k = 1;
            end
            for i = ii % connect with all other nodes
                for j = jj
                    if (k == 0) || (i ~= j)               %(k == 0) || (i ~= j)  %(abs(i) ~= abs(j))
                        rx = i; %0 oder 1
                        ry = j; % -1, 0 , 1
                        if (1 <= elx+rx) && (elx+rx <= nelx ) && (1 <= ely+ry) && (ely+ry <= nely )
                            k2 = (ely+ry)+nely*((elx+rx)-1); % node number
                            if (k1 < k2) && (( sn == 1) || (( sn > 1) && (( staebe(sn-1,1) ~= k1) || ( staebe(sn-1 ,2) ~= k2))))
                                staebe(sn,1) = k1; 
                                staebe(sn,2) = k2; 
                                sn = sn+1;
                            end
                        end
                    end
                end % for j
            end % for i
        end % for ely
    end % for elx
end % for elz  
    
elseif art == 3 || art ==9 || art==10 % 6 Ecke
    
    for elx = 1:nelx
        for ely = 1:nely
            k1 = ely+nely*(elx-1) ; % node number
            %Waagrechte Stäbe 1
            if mod(elx-2,4) == 0
                if mod(ely,2) == 1
                    if elx ~= nelx
                        k2 = k1 + nely;
                        staebe(sn,1) = k1;
                        staebe(sn,2) = k2;
                        sn = sn+1;
                    end
                end
            end
            %Waagrechte Stäbe 2
            if mod(elx,4) == 0
                if mod(ely,2) == 0
                    if elx ~= nelx
                        k2 = k1 + nely;
                        staebe(sn,1) = k1;
                        staebe(sn,2) = k2;
                        sn = sn+1;
                    end
                end
            end
            %Schräge Stäbe mit Steigung 1
            if mod(elx-1,4) == 0
                if mod(ely,2) == 0
                    if elx ~= nelx && ely ~= nely
                        k2 = k1 + nely +1;
                        staebe(sn,1) = k1;
                        staebe(sn,2) = k2;
                        sn = sn+1;
                    end
                end
            end
            %Schräge Stäbe mit Steigung 1
            if mod(elx+1,4) == 0
                if mod(ely-1,2) == 0
                    if elx ~= nelx && ely ~= nely
                        k2 = k1 + nely +1;
                        staebe(sn,1) = k1;
                        staebe(sn,2) = k2;
                        sn = sn+1;
                    end
                end
            end
            %Schräge Stäbe mit Steigung -1
            if mod(elx,4) == 0
                if mod(ely,2) == 0
                    if  ely ~= nely
                        k2 = k1 - nely +1;
                        staebe(sn,1) = k1;
                        staebe(sn,2) = k2;
                        sn = sn+1;
                    end
                end
            end
            %Schräge Stäbe mit Steigung -1
            if mod(elx+2,4) == 0
                if mod(ely,2) == 1
                    if  ely ~= nely
                        k2 = k1 - nely +1;
                        staebe(sn,1) = k1;
                        staebe(sn,2) = k2;
                        sn = sn+1;
                    end
                end
            end
        end % for ely
    end % for elx
    
end    


%Zusätzliche Streben in Wabenstruktur für Gitter 10,11 
rtol = 1e-4; %Relativ Toleranz für Knoten Erkennung

if art==10
    for k=3:2:(nely-1)
        staebe(sn,1)=k;
        staebe(sn,2)=k+nely;
        sn=sn+1;
        if nely/2 == round(nely/2)
            staebe(sn,1)=nely*nelx-k+1;
            staebe(sn,2)=nely*nelx-k-nely+1;
            sn=sn+1;
        else
            staebe(sn,1)=nely*nelx-k;
            staebe(sn,2)=nely*nelx-k-nely;
            sn=sn+1;
        end
    end
    if nely/2 ~= round(nely/2)
        staebe(sn,1)=nely*nelx-1;
        staebe(sn,2)=nely*nelx-nely-1;
        sn=sn+1;
    end
end

if art==9 || art==10
    for k =(nely*nelx+1) : length(knoten)
        if art==10 || knoten(k,2) == gy0 || knoten(k,2) == gy1
            if ismembertol([knoten(k,1)-1.5*((gx1-gx0)/(nelx-1)), knoten(k,2)], knoten, 'ByRows',rtol*((gx1-gx0)/(nelx-1))) % Strebe nach links
                [~,k2]=ismembertol([knoten(k,1)-1.5*((gx1-gx0)/(nelx-1)), knoten(k,2)], knoten, 'ByRows', rtol*((gx1-gx0)/(nelx-1)));
                staebe(sn,1) = k;
                staebe(sn,2) = k2;
                sn = sn+1;
            end
            
            if ismembertol([knoten(k,1)+1.5*((gx1-gx0)/(nelx-1)), knoten(k,2)], knoten, 'ByRows', rtol*((gx1-gx0)/(nelx-1))) % Strebe nach rechts
                [~,k2]=ismembertol([knoten(k,1)+1.5*((gx1-gx0)/(nelx-1)), knoten(k,2)], knoten, 'ByRows', rtol *((gx1-gx0)/(nelx-1)));
                staebe(sn,1) = k;
                staebe(sn,2) = k2;
                sn = sn+1;
            end
        end
        if ismembertol([knoten(k,1)-0.5*((gx1-gx0)/(nelx-1)), knoten(k,2) + ((gy1-gy0)/(nely-1))], knoten, 'ByRows', rtol*min([((gx1-gx0)/(nelx-1)),((gy1-gy0)/(nely-1))])) % Strebe nach links oben
            [~,k2]=ismembertol([knoten(k,1)-0.5*((gx1-gx0)/(nelx-1)), knoten(k,2)+((gy1-gy0)/(nely-1))], knoten, 'ByRows', rtol*min([((gx1-gx0)/(nelx-1)),((gy1-gy0)/(nely-1))]));
            staebe(sn,1) = k;
            staebe(sn,2) = k2;
            sn = sn+1;
        end
        
        if ismembertol([knoten(k,1)+0.5*((gx1-gx0)/(nelx-1)), knoten(k,2) + ((gy1-gy0)/(nely-1))], knoten, 'ByRows', rtol*min([((gx1-gx0)/(nelx-1)),((gy1-gy0)/(nely-1))])) % Strebe nach rechts oben
            [~,k2]=ismembertol([knoten(k,1)+0.5*((gx1-gx0)/(nelx-1)), knoten(k,2)+((gy1-gy0)/(nely-1))], knoten, 'ByRows', rtol*min([((gx1-gx0)/(nelx-1)),((gy1-gy0)/(nely-1))]));
            staebe(sn,1) = k;
            staebe(sn,2) = k2;
            sn = sn+1;
        end
        
        if ismembertol([knoten(k,1)-0.5*((gx1-gx0)/(nelx-1)), knoten(k,2) - ((gy1-gy0)/(nely-1))], knoten, 'ByRows', rtol*min([((gx1-gx0)/(nelx-1)),((gy1-gy0)/(nely-1))])) % Strebe nach links unten
            [~,k2]=ismembertol([knoten(k,1)-0.5*((gx1-gx0)/(nelx-1)), knoten(k,2)-((gy1-gy0)/(nely-1))], knoten, 'ByRows', rtol*min([((gx1-gx0)/(nelx-1)),((gy1-gy0)/(nely-1))]));
            staebe(sn,1) = k;
            staebe(sn,2) = k2;
            sn = sn+1;
        end
        
        if ismembertol([knoten(k,1)+0.5*((gx1-gx0)/(nelx-1)), knoten(k,2) - ((gy1-gy0)/(nely-1))], knoten, 'ByRows', rtol*min([((gx1-gx0)/(nelx-1)),((gy1-gy0)/(nely-1))])) % Strebe nach rechts unten
            [~,k2]=ismembertol([knoten(k,1)+0.5*((gx1-gx0)/(nelx-1)), knoten(k,2)-((gy1-gy0)/(nely-1))], knoten, 'ByRows', rtol*min([((gx1-gx0)/(nelx-1)),((gy1-gy0)/(nely-1))]));
            staebe(sn,1) = k;
            staebe(sn,2) = k2;
            sn = sn+1;
        end
    end
end




%Ränder für Wabenstrukturen
if art == 3 || art==7 || art==8 || art==9 || art==10
    for elx = 1:nelx
        for ely = 1:nely
            k1 = ely+nely*(elx-1) ; % node number
            
            if elx == 1 || elx == nelx %linker und rechter Rand
                if ely ~= nely
                    k2 = k1 + 1;
                    if ~ismember(k2, staebe) && k2~=nely && k2 ~=nely*nelx
                        k2 = k1 + 2;
                    end
                    if (ismember(k1, staebe) || k1 == 1 || k1==(nelx-1)*nely+1)
                        staebe(sn,1) = k1;
                        staebe(sn,2) = k2;
                        sn = sn+1;
                    end
                end
            end
            
            if ely == 1 || ely == nely %oberer und unterer Rand
                if elx ~= nelx
                    if ismember(k1, staebe)
                        k2 = k1 + nely;
                        if ~ismember(k2, staebe) && k2+nely<= nely*nelx
                            k2 = k1 + 2*nely;
                            if ~ismember(k2, staebe) && k2+nely<= nely*nelx
                                if art==9 || art==10
                                    continue;
                                else
                                    k2 = k1 + 3 * nely;
                                end
                            end
                        end
                        staebe(sn,1) = k1;
                        staebe(sn,2) = k2;
                        sn = sn+1;
                    end
                end
            end
            
        end % for ely
    end % for elx
end


staebe = staebe(1:sn-1,:);
size( staebe );
staebe = unique(staebe,'rows');
size( staebe );




if art==3 || art==7 || art==8  || art==9 || art==10 %Löschen ungenutzer Knoten bei 6-Ecken -> Mappen auf neue Knotenmatrix
    u_kn=unique([staebe(:,1); staebe(:,2)]);
    u_kn=u_kn(1:end,:);
    u_kn_koord = knoten(u_kn,:);
    [~,index] = ismember(u_kn_koord,knoten,'rows');
    knoten=u_kn_koord;
    M = containers.Map(index,[1:length(u_kn)]);
    for i=1 : length(staebe)
        for l=1 :2
            staebe(i,l) = M(staebe(i,l));
            
        end
    end
    
end



% Aufbringen der zufälligen Knotenverschiebung für 6-Eck Strukturen
if ismember(art,[3,7,8,9,10])
    for i=1:length(knoten)
        if knoten(i,1)~=gx0 && knoten(i,2)~=gy0 && knoten(i,1)~=gx1 && knoten(i,2)~=gy1
            knoten(i,1) = knoten(i,1) + varx * randi([-1,1]) * rand(1);
            knoten(i,2) = knoten(i,2) + vary * randi([-1,1]) * rand(1);
        end
    end
end





end

function [z,n]= kuerze(z,n)
% Die Zahlen z und n werden durch den ggT geteilt .
if (n < 0) && (z < 0)
    z=-z;
    n=-n;
end
if z == 0
    n = 1;
elseif n == 0
    z = 1;
elseif (z ~= 1) && (n ~= 1)
    %ggt=gcd(z,n);
    ggt=ggT(z,n);
    z = z/ggt;
    n = n/ggt;
end
end

function ggt=ggT(z,n)
% Der groesste gemeinsame Teiler ggT von z und n wird mittels Euklid'schenAlgorithmus bestimmt.
% Die Funktion gcd von MatLab ist langsamer und maechtiger .
z = abs(z);
n = abs(n);
while (n ~= 0)
    tmp = n;
    n = z - floor (z/n) * n;
    z = tmp;
end
ggt = z;
end

