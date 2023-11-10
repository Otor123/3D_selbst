%function trab_evo
function[e_mat,A,dof,V,P,S_mac,sigma]=...
    trab_evo(q0,edof,u_pre,mat,nel,ndof,a_alt,F_mac,e_mat)
%MESO FE calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
tol    = 1e-7;                   % tolerance of newton iteration
%e_mat  = extr_dof(edof,q0);       % init material coordinates
dof    = q0;                      % init spatial coordinates       
a_var  = a_alt;                      % init internal variables %Querschnitt
F_pre  = zeros(ndof,1);             % init incremental neuman bc's
u_inc  = u_pre;                   % init incremental dirichlet bc's
lodfac = 0;                       %init loadfactor increments with every nlod(is)=1 entry
for ie=1:nel emat(ie) = 1; end;
nip=1;
iter=0;

[n,nedof] = size(edof);
nedof     = nedof-1;
residuum = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% global newton-raphson iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    while residuum > tol
      iter=iter+1;

      A = zeros(nel,nip);             % init internal variables
      R = zeros(ndof,1);              % init global residual
      K = sparse(ndof,ndof);          % init global stema
      e_spa = extr_dof(edof,dof);     % extract dofs of global vector
      e_par = mat(emat(:),:);         % extract material parameters
      V=0;

	  kkv = zeros(nedof*nedof,nel);   % matrix coefficients
	  kki = zeros(nedof*nedof,nel);   % matrix i indexes
	  kkj = zeros(nedof*nedof,nel);   % matrix j indexes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%% loop over all elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nel aus rve_ini.m -> nel = size(dof,1) (degrees of freedome of every node)

      for ie = 1:nel %anzahl der stäbe
            if nedof==4
	    [Ke,Re,Ae,Ve] = stab2d2(e_mat(ie,:),e_spa(ie,:),a_var(ie),e_par(ie,:));                
        elseif nedof==6 
	    [Ke,Re,Ae,Ve] = stab3d3(e_mat(ie,:),e_spa(ie,:),a_var(ie,:),e_par(ie,:));
		end
%%%%%%%% assembly with sparse solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%     	[K, R, I ] = assm_sys(edof(ie,:),K,Ke,R,Re,I,Ie);
		index       = edof(ie,2:end);
		rows		= index(1:nedof)'*ones(1,nedof);
		cols		= ones(nedof,1)*index(1:nedof);
		kkv(:,ie) = Ke(:);
		kki(:,ie) = rows(:);
		kkj(:,ie) = cols(:);
		R(edof(ie,2:end)) = R(edof(ie,2:end)) + Re;
		A(ie,:)=Ae(:);
        %x_tna(ie,:)=x_tn;
        V=V+Ve;
	  end
	  K = sparse(kki(:),kkj(:),kkv(:),ndof,ndof);  % gesamt Steifigkeit stäbe  
%%%%% loop over all elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%wir haben nur ein ellement, deswegen keine schleife%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      
      F_inc      = lodfac * F_pre;        % scale neumann   bc's with dt
      if mean(F_pre)==0
          u_inc(:,2)=u_pre(:,2);
      end
        
      R = R - F_inc;                    % add ext load to righthand side
      dofold = dof; 
      [dof,F] = solve_nr(K,R,dof,iter,u_inc);     % solution and  update 
      residuum=res_norm((dof-dofold),u_inc);

      if length(dof)==4 || length(dof)==6
          if iter>1
            residuum=0;
          else
              residuum=1;
          end
      end

  %    disp(['*** meso iteration: ',num2str(iter),...
   %         ' *** meso residual: ',num2str(residuum)]);  
	  if(iter>300); 
          disp(['*** NO GLOBAL MESO CONVERGENCE ***']);
          residuum
              [P,sigma]=stress_hom(F,q0,dof,u_pre);
                S_mes=stiff_red(K,u_pre,ndof);
                S_mac=stiff_mac(S_mes,dof,q0,R,sigma,F_mac,u_pre);
          return; 
      end;		    
end

    %Berechnung der homogenisierten Spannung und Steifigkeit
         if nedof==4
	    [P,sigma]=stress_hom(F,q0,dof,u_pre);             
        elseif nedof==6 
	    [P,sigma]=stress_hom3D(F,q0,dof,u_pre);
		end

    
    S_mes=stiff_red(K,u_pre,ndof);
    
        if nedof==4
	    S_mac=stiff_mac(S_mes,dof,q0,R,sigma,F_mac,u_pre);              
        elseif nedof==6 
	    S_mac=stiff_mac3D(S_mes,dof,q0,R,sigma,F_mac,u_pre);  
		end
%%% global newton-raphson iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
