%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q0,edof,emat,bc,F_ext,mat,ndim,nel,node,ndof,nip,nlod] = ex_square_test(F); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% material parameters for density growth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emod = 10000;             %sollte der gleiche wie im mesomodell sein
nue  = 0.33;      %Poisonratio?
rho0 =0.4;        %wird in mesomodell überführt
psi0 = 0.0014;    %ist laut jacob's studien im schnitt etwa 7 mal kleiner als im stabmodell
expm = 3.0;       %-
expn = 2.0;       %-
dt   = 0.01;       %entscheidend für konvergenz
c_rs = 200;       %ist laut jacob's studien im schnitt etwa 20 mal größer als im stabmodell

mat = [emod,nue,rho0,psi0,expm,expn,dt,c_rs];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate node points q0 and connectivity edof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xbox(1) = 0.0;
xbox(2) = 1;
ybox(1) = 0.0;
ybox(2) = 1;                
zbox(1) = 0.0;
zbox(2) = 1;   

nx =  1; %Teilung auf X-Achse
ny =  1; %Teilung auf Y-Achse
nz =  1; %Teilung auf Y-Achse

[q0,edof] = mesh_cub(xbox,ybox,zbox,nx,ny,nz);  

[nel, sizen] = size(edof);
[ndof,sizen] = size(q0);
 for ie=1:nel emat(ie) = 1; end;
 node        = ndof/2;
 nip  = 2; %4
 ndim = 3; %2
 
if F==1 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dirichlet boundary conditions // sind immer Lagerungen
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bc(1,1)= 1;   bc(1,2)=0;                  
    bc(2,1)= 2;   bc(2,2)=0;
    bc(3,1)= 3;   bc(3,2)=0;
    bc(4,1)= 4;   bc(4,2)=0;                  
    bc(5,1)= 5;   bc(5,2)=0;
    bc(6,1)= 6;   bc(6,2)=0;
    bc(7,1)= 7;   bc(7,2)=0;
    bc(8,1)= 8;   bc(8,2)=0;
    bc(9,1)= 9;  bc(9,2)=0;
    bc(10,1)= 10; bc(10,2)=0;
    bc(11,1)= 11; bc(11,2)=0;
    bc(12,1)= 12; bc(12,2)=0;
    bc(13,1)= 13; bc(13,2)=0;
    bc(14,1)= 15; bc(14,2)=0;
    bc(15,1)= 16; bc(15,2)=0;
    bc(16,1)= 18; bc(16,2)=0;
    bc(17,1)= 19; bc(17,2)=0;
    bc(18,1)= 21; bc(18,2)=0;
    bc(19,1)= 22; bc(19,2)=0;
    bc(20,1)= 24; bc(20,2)=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neumann boundary conditions // sind immer die Kräfte
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F_ext =  zeros(ndof,1);
    
%     F_ext(4) = -0.5; %y-direction k2
%     F_ext(8) = -0.5; %y-direction k4

    F_ext(14) = -0.5; %y-direction k2
    F_ext(17) = -0.5; %y-direction k4
    F_ext(20) = -0.5; %y-direction k6
    F_ext(23) = -0.5; %y-direction k8   
    
    %last wird stufenweise
    d_ges=8;
    nlod=zeros(d_ges/dt,1);
    nlod(1)=1;
    nlod(2/mat(7))=1;
    nlod(4/mat(7))=1;
    nlod(6/mat(7))=1;
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dirichlet boundary conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     bc(1,1)= 2; bc(1,2)=0;                  
%     bc(2,1)= 4; bc(2,2)=0;
%     bc(3,1)= 5; bc(3,2)=0;
%     bc(4,1)= 6; bc(4,2)=0;                  
%     bc(5,1)= 7; bc(5,2)=0;
%     bc(6,1)= 8; bc(6,2)=0;
    bc(1,1)= 2;     bc(1,2)=0;
    bc(2,1)= 3;     bc(2,2)=0;     
    bc(3,1)= 5;     bc(3,2)=0;
    bc(4,1)= 6;     bc(4,2)=0;     
    bc(5,1)= 7;     bc(5,2)=0;
    bc(6,1)= 8;     bc(6,2)=0; 
    bc(7,1)= 9;     bc(7,2)=0;     
    bc(8,1)= 10;    bc(8,2)=0;
    bc(9,1)= 11;    bc(9,2)=0;
    bc(10,1)= 12;   bc(10,2)=0;     
    bc(11,1)= 14;   bc(11,2)=0;  
    bc(12,1)= 15;   bc(12,2)=0;     
    bc(13,1)= 17;   bc(13,2)=0;
    bc(14,1)= 18;   bc(14,2)=0;     
    bc(15,1)= 19;   bc(15,2)=0;
    bc(16,1)= 20;   bc(16,2)=0;  
    bc(17,1)= 21;   bc(17,2)=0;     
    bc(18,1)= 22;   bc(18,2)=0;
    bc(19,1)= 23;   bc(19,2)=0; 
    bc(20,1)= 24;   bc(20,2)=0; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neumann boundary conditions%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %In jedem Punkt kann eine Kraft angreifen!! (Bei nx=2 und ny=2 -> 18
    %Einträge // bei nx = 1 und ny = 1 -> 8 Einträge
    F_ext =  zeros(ndof,1);

%     F_ext(1)  = 0.5; %x-direction k1
%     F_ext(3)  = 0.5; %x-direction k2
    
    F_ext(1)  = 0.5; %x-direction k1
    F_ext(4)  = 0.5; %x-direction k2
    F_ext(13) = 0.5; %x-direction k5
    F_ext(16) = 0.5; %x-direction k6    

    %last wird stufenweise (aktuell bei 1, 200, 400, 600)
    d_ges=8;
    nlod=zeros(d_ges/dt,1);
    nlod(1)=1;
    nlod(2/mat(7))=0.5;
    nlod(4/mat(7))=1;
    nlod(6/mat(7))=1;    
%     nlod(1)=2;
%     nlod(2/mat(7))=1.5;
%     nlod(4/mat(7))=2;
%     nlod(6/mat(7))=2;
end
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nmat,ngrw] = size(mat);
if (ngrw==8); disp('*** density growth ***'); end;  % density
if (ngrw==9); disp('*** volume growth ***' ); end;  % volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%