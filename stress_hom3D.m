function [P,sig]=stress_hom3D(F,q0_k,dof_k,u_pre)
%2D-Berechnung der homogenisierten Spannung
%finden der Randknoten
rand1=u_pre(:,1);
dof=dof_k(rand1);
q0=q0_k(rand1);
R=F(rand1);
rand=u_pre(2:2:end,1)/2;
nr=length(rand);
%Aufteilung in x und y und z Koordinaten
dofx=dof(1:3:length(q0));
dofy=dof(2:3:length(q0));
dofz=dof(3:3:length(q0));
q0x=q0(1:3:length(q0));
q0y=q0(2:3:length(q0));
q0z=q0(3:3:length(q0));
o=find(q0y==(-q0y(1)));
u=find(q0y==(q0y(1)));
lo=o(find(q0x(o)==q0x(1)));
ro=o(find(q0x(o)==-q0x(1)));
ru=u(find(q0x(o)==-q0x(1)));


%Aufteilung der Reaktionskräfte in x und y und z
F_reac(:,1)=R(1:3:length(R));
F_reac(:,2)=R(2:3:length(R));
F_reac(:,3)=R(3:3:length(R));

%Berechnung der RVE-Fläche
V_sq=polyarea([q0x(1) q0x(lo) q0x(ro) q0x(ru)],[q0y(1) q0y(lo) q0y(ro) q0y(ru)]);

V_def=polyarea([dofx(1) dofx(lo) dofx(ro) dofx(ru)],[dofy(1) dofy(lo) dofy(ro) dofy(ru)]);
if length(q0)==6
    V_sq=1;
    V_def=1*(1-2*u_pre(2,2));
end

%Randkoordinaten
randknot(:,1)=q0x;
randknot(:,2)=q0y;
randknot(:,3)=q0z;

randknotd(:,1)=dofx;
randknotd(:,2)=dofy;
randknotd(:,3)=dofz;
[randis,dim]=size(randknot);

%Berechnung der homogenisierten Kirchhoff Spannung
Pe=0;
for i=1:randis
    Pe=Pe+F_reac(i,:)'*randknot(i,:);
end
P=Pe/V_sq;

%Berechnung der Cauchy Spannung
sige=0;
for i=1:randis
    sige=sige+F_reac(i,:)'*randknotd(i,:);
end
sig=sige/V_def;








