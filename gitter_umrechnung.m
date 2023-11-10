function [u_pre] = gitter_umrechnung(q0,F,randis)
%Umrechnung zurück von q0 und edof zu knoten und Stäbe
%dann Umrechnung von F auf u_pre
for i=1:(size(q0,1)/2)
    knoten(i,1)=q0(2*i-1,1); %x-werte
    knoten(i,2)=q0(2*i,1); %y-werte
end


randknoten=knoten(randis,:);     
z=length(randis);

for m=1:z
    versch=F*randknoten(m,:)'-randknoten(m,:)';     %(Berechnung Verschiebung der Randknoten)Verschiebung=Deformationsgradient F multipliziert mit den undef. Koordinaten- undef. Koordinaten
    %in x-Richtung
    ver(2*m-1,1)=randis(m)*2-1;     %Freiheitsgrad der Verschiebung nur der Randknoten
    ver(2*m-1,2)=versch(1);         %Wert der Verschiebung
    %in y-Richtung
    ver(2*m,1)=randis(m)*2;
    ver(2*m,2)=versch(2);
end


[a,b]=sort(ver(:,1));
u_pre=ver(b,:);
u_pre(:,2)=u_pre(:,2);
end

