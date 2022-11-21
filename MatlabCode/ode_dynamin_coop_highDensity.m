 % Dsol0 is in units of uM
 % Dmem0, Rclus0 and Rdilute0 are in units of /um2!
% Length scale h is in units of um (e.g. 0.001um)
% time is in s  it is either [T0 Tfinal] or an array of time pts.  
% 
% k in units of uM(-1)s(-1)
%kb in s(-1)
%VAratio is units of um
%Add in a reaction where dynamin within the cluster recruits dynamin. 
function [timepts, conc] = ode_dynamin_coop_highDensity(V, Atotal, Acluster,h, time,...
    densityIncrease,...
    kfDRC, Rdilute0, Dsol0_1, Dmem0, ...
    kfDRD, kfDydy, Rclus0) 


DynTotal=Dsol0_1+Dmem0*Atotal/V/602.0; %total Dynamin in solution and on the membrane
%in previous simulation
Dmem0=Dmem0*densityIncrease; %the copies on the membrane are now higher
Dmem0=Dmem0*Atotal/V/602; %now in uM
%cluster density and recruiter (unbound) density is considered unchanged.
Rclus0=Rclus0*Acluster/V/602; %now in uM
Rdilute0=Rdilute0*Atotal/V/602; %now in uM
VAratio=V/Atotal;
display('Starting solution conc of Dynamin')
Dsol0=DynTotal-Dmem0
kbDRD=Dsol0*Rdilute0/Dmem0*kfDRD;%initially, the dynamin is in equilibrium with the membrane
%make all cluster binding irreversible
kbDRC=0; %unbinding from the cluster 
KDdilute=kbDRD/kfDRD;
display('off rate')
kbDRD
display('KD (uM)')
KDdilute
display('gamma ')
V/Atotal/h

display('rate to bind cluster')
kfDRC

y0=zeros(6,1);
y0(1)=Dsol0;%dynamin in solution 
y0(2)=Dmem0;%dynamin on membrane
y0(3)=Rclus0; %recruiter on the clsuter membrane
y0(4)=Rdilute0;%recruiter on the membrane dilute
y0(5)=0;%DmemR, dynamin from 2D +recruiter in cluster
y0(6)=0; %DsolR, dynamin from 3D+recruiter in cluster

gamma=VAratio/h;% dimensionless!

opt=odeset('RelTol',1E-6,'AbsTol',1E-7);
[timepts,conc] = ode23s(@(t,y) odes_6(t,y, kfDRC, kfDRD,  kbDRC, kbDRD, gamma, kfDydy),time,y0, opt);




function dy = odes_6(t,y, kfDRC, kfDRD,  kbDRC, kbDRD, gamma, kfDydy)

dy = zeros(6,1);    % a column vector

bind3DtoCluster=-kfDRC*y(1)*y(3);
bind2DtoCluster=-kfDRC*gamma*y(2)*y(3);
bind3DtoDilute=-kfDRD*y(1)*y(4);
bind2DtoDyn = -kfDydy*gamma*y(2)*y(5);
bind2Dto3DDyn = -kfDydy*gamma*y(2)*y(6);
bind3DtoDyn = -kfDydy*y(1)*y(5);
bind3Dto3DDyn = -kfDydy*y(1)*y(6);
%DYN-DYN interactions, only can occur within the cluster
%y(2)+y(5)->2y(5)  ; 2D
%y(2)+y(6)->y(5)+y(6) ; 2D
%y(1)+y(5)->y(5)+y(6) ; 3D
%y(1)+y(6)->2y(6) ; 3D

%dynamin in solution
dy(1) = bind3DtoCluster+kbDRC*y(6)+bind3DtoDilute+kbDRD*y(2)+bind3DtoDyn+bind3Dto3DDyn;
%dynamin on membrane
dy(2) = bind2DtoCluster+kbDRC*y(5)-bind3DtoDilute-kbDRD*y(2)+bind2DtoDyn+bind2Dto3DDyn;
%recruiter on the cluster membrane
dy(3) = bind3DtoCluster + bind2DtoCluster +kbDRC*y(5)+kbDRC*y(6);  
%recruiter on the membrane dilute
dy(4) = bind3DtoDilute+kbDRD*y(2); 
%Dynamin in cluster from 2D
dy(5) = -bind2DtoCluster-kbDRC*y(5)-bind2DtoDyn-bind2Dto3DDyn;
%Dynamin in cluster from 3D
dy(6) = -bind3DtoCluster - kbDRC*y(6)-bind3DtoDyn-bind3Dto3DDyn;




