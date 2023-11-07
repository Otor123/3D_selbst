%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [var,dres]=updt_2a(epsi,var,mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% density update for functional adaption
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:  F        -> deformation gradient
%         var      -> var = a / a_0 - 1 
%         mat      -> material parameters
% output: var      -> var = a / a_0 - 1 
%         facs     -> factor for stresses    
%         fact     -> factor for tangent operator    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-12; 

emod = mat(1);   a0 = mat(2);   Wa = mat(3);   
expm = mat(4);   expn = mat(5);   dt   = mat(6); c_rs=mat(7);


W=0.5*emod*epsi^2; %Formel 5.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% euler backward - implicit time integration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
a_k0 = var;
a_k1 = var; 

if a_k0 > 1e-12 %wenn stab noch nicht aufgelöst ist
%if W > 1e-12
iter = 0;   res  = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% local newton-raphson iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    while abs(res) > tol
      iter=iter+1;
  	  
      res =a_k1-a_k0-c_rs*a_k1*(W-Wa)*dt ;%new W-dehnungsenergie Wa-vorgegebener Wert 
	  dres=1-c_rs*(W-Wa)*dt; %Ableitung von res
       
      da  =-res/dres ; %new area
	  a_k1 = a_k1+da;
 
	  if(iter>80) 
          disp(['*** NO LOCAL MESO CONVERGENCE ***']); 
          return; 
      end
    end
%%% local newton-raphson iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if a_k1<1e-12 %Wenn stab quasi aufgelöst ist
    a_k1=1e-12;
end
else %stab ist aufgelöst
    a_k1=a_k0;
end

var = a_k1;

dres= 1-c_rs*(W-Wa)*dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%