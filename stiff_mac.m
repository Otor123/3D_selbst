function [S_mac]=stiff_mac(S_mes,dof_k,q0_k,R_k,sigma,F_mac,u_pre)
%Berechnung der macroskopischen Steifigkeitsmatrix
%alle dofs der Randknoten müssten der ersten Spalte von u_pre entsprechen
S_mes=full(S_mes);
rand1=u_pre(:,1);
dof=dof_k(rand1);
q0=q0_k(rand1);
R=R_k(rand1);
rand=u_pre(2:2:end,1)/2;
nr=length(rand);

% dofx=dof(1:2:length(q0));
% dofy=dof(2:2:length(q0));
% %Fläche A
% A_def=(max(dofx)-min(dofx))*(max(dofy)-min(dofy));

dofx=dof(1:2:length(q0));
dofy=dof(2:2:length(q0));
q0x=q0(1:2:length(q0));
q0y=q0(2:2:length(q0));
o=find(q0y==(-q0y(1)));
u=find(q0y==(q0y(1)));
lo=o(find(q0x(o)==q0x(1)));
ro=o(find(q0x(o)==-q0x(1)));
ru=u(find(q0x(o)==-q0x(1)));
%Fläche A
%A_def=polyarea([dofx(1) dofx(lo) dofx(ro) dofx(ru)],[dofy(1) dofy(lo) dofy(ro) dofy(ru)]);
A_def=polyarea([q0x(1) q0x(lo) q0x(ro) q0x(ru)],[q0y(1) q0y(lo) q0y(ro) q0y(ru)]);
if length(q0)==4
    A_def=1;
end

S_m=0;
for i=1:nr
    su=0;
    for j=1:nr
        su=su+outer(S_mes(i*2-1:i*2,j*2-1:j*2),q0(j*2-1:j*2));
    end
    S_m=S_m+outer(q0(i*2-1:i*2),su);
end
S_m2=permute(S_m,[2 1 3 4]);
S_mac=1/A_def*S_m2;

% S_m=0;
% for i=1:nr
%     su=0;
%     for j=1:nr
%         su=su+outer(S_mes(i*2-1:i*2,j*2-1:j*2),q0(j*2-1:j*2));
%     end
%     S_m=S_m+outer(dof(i*2-1:i*2),su)+outer(R(i*2-1:i*2),outer(eye(2),q0(i*2-1:i*2)));
% end
% S_mac=1/A_def*S_m;




