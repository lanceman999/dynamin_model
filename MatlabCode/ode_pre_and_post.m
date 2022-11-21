%Solve for the conc/copies of dynamin pre and post stimulation.
%Read in a vector of 8 parameters, that were solved for using the Python
%genetic algorithms:
%read in optimal parameters, assign to values
%INPUTS: (these are the first 8 values printed in optFitParms_prePost.txt
%files).
% kfDRC=parmVector8(1);
% Rdilute0=parmVector8(2);
% Dsol0=parmVector8(3);
% Dmem0=parmVector8(4);
% kfDRD=parmVector8(5);
% kfDydy=parmVector8(6); 
% Rclus0=parmVector8(7);
% kfDense=parmVector8(8);

%Input a string for isoform: either 'BB','AB', or 'AA'
%OUTPUTS:
%time points for solutions (Nx1)
%concPre: (Nx6), concentrations in uM for all species pre stimulation
%concP0st: (Nx6), concentrations in uM for all species post stimulation
%conc(:,1)=Dsol
%conc(:,2)=Dmem
%conc(:,3)=Aclus
%conc(:,4)=Rdilute
%conc(:,5)=Dclus2D
%conc(:,6)=Dclus3D 
%copies=conc*V*602
%total dynamin in cluster is dynamin_copies/SA *Acluster:
%Acluster*(Dclus2D_copy/Acluster+Dclus3D_copy/Acluster+Dmem_copy/Atotal);
% where we note that Dmem is on the whole surface, including in the cluster

function[timepts, concPre, concPost]=ode_pre_and_post(parmVector8,isoform)

%time points range, and to print solution
time=logspace(-5,0.65,500);
h=0.01; %um. 3D to 2D equilibrium constant change
VAratio=1.9 ; %um 

if ~exist('isoform')
    isoform='BB'
end

if(isoform == 'BB')
    %%%%%%%%%%%%%%%%%%%Parameters for Dyn1BB
    V=3.53; %um3
    Apre=1.85 ; %um2
    Aclus=0.022 ; %um2
    Vpost=0.615; %um3
    Apost=0.322; %um2
    AclusPost=0.028; %um2
    densityIncrease = 5;
    targetPre = 30;
    targetPost=10;
    display(' Simulating BB !')
elseif(isoform == 'AB')
        %%%%%%%%%%%%%%%%%%%Parameters for Dyn1AB
         V=8.25; %um3
         Apre=4.34 ; %um2
         Aclus=0.023 ; %um2
         Vpost=1.083; %um3
         Apost=0.57; %um2
         AclusPost=0.026; %um2
         densityIncrease = 7; %factor increase
        targetPre = 45;
        targetPost=14;
        display(' Simulating AB !')
elseif(isoform =='AA')      
    %%%%%%%%%%%%%%%%%%%Parameters for Dyn1AA
    Apre=2.85 ; %um2
    V=Apre*VAratio; %um3
    Aclus=0.02 ; %um2
    Apost=1.43; %um2
    Vpost=Apost*VAratio; %um3
    AclusPost=0.02; %um2
    densityIncrease = 2;
    targetPre = 20;
    targetPost=23;
    display(' Simulating AA !')
end

%read in optimal parameters, assign to values
kfDRC=parmVector8(1);
Rdilute0=parmVector8(2);
Dsol0=parmVector8(3);
Dmem0=parmVector8(4);
kfDRD=parmVector8(5);
kfDydy=parmVector8(6); 
Rclus0=parmVector8(7);
kfDense=parmVector8(8);

%Solve dynamin pre-stimulation
[timepts, concPre] = ode_dynamin_coop_density(V, Apre, Aclus,h, time,...
    kfDRC, Rdilute0, Dsol0, Dmem0, ...
    kfDRD, kfDydy, Rclus0);

copies=concPre*602*V;%conc is in units of uM. convert to copies.
acd=copies(:,5)/Aclus;
ad=copies(:,6)/Aclus;
display('pre Stim 2D/3D recruitment (final)')
acd(end)/ad(end)
display('pre Stim 2D/3D recruitment (initial)')
acd(2)/ad(2)
display('pre Stim 2D/3D recruitment (mean)')
0.5*acd(end)/ad(end)+0.5*acd(2)/ad(2)

acPre=concPre(:,2)*602*V/Apre*Aclus; %copies in the cluster
%solve dynamin post-stimulation. Initial density is higher, and rate of
%binding activator can change. Will correct for decrease in solution copies
%of dynamin

display('POST STIMULATION RESULTS')
[timepts, concPost] = ode_dynamin_coop_highDensity(Vpost, Apost, AclusPost,h, time,...
    densityIncrease,...
    kfDense, Rdilute0, Dsol0, Dmem0, ...
    kfDRD, kfDydy, Rclus0);

copies=concPost*602*Vpost;%conc is in units of uM. convert to copies.
acd=copies(:,5)/AclusPost;
ad=copies(:,6)/AclusPost;
display('POST Stim 2D/3D recruitment (final)')
acd(end)/ad(end)
display('POST Stim 2D/3D recruitment (initial)')
acd(2)/ad(2)
display('POST Stim 2D/3D recruitment (mean)')

0.5*acd(end)/ad(end)+0.5*acd(2)/ad(2)
acPost=concPost(:,2)*602*Vpost/Apost*AclusPost; %copies in the cluster
edat=load('expdataBB_from4.dat');

%plot pre
relPre=plot_dynamin_fromODE(time, concPre, V, Apre, Aclus, 1);
del1=(edat(end,2)-targetPre*edat(1,2))/(1.0-targetPre);
display('del Pre: ')
del1

figure(3)
plot(edat(:,1)+4,(edat(:,2)-del1)/(edat(1,2)-del1)*acPre(1),'r-', 'LineWidth',2);
%plot post
del2=(edat(end,2)-targetPost*edat(1,2))/(1.0-targetPost);
display('del Post: ')
del2
relPost=plot_dynamin_fromODE(time, concPost, Vpost, Apost, AclusPost, 4);
figure(6)
plot(edat(:,1)+4,(edat(:,2)-del2)/(edat(1,2)-del2)*acPost(1),'c-','LineWidth',2);

%PLOT BOTH SOLUTIONS (PRE AND POST) ON ONE AXIS
f3=figure(7);
 ax3=axes('Parent',f3,'FontSize',20,'LineWidth',1);
hold(ax3);
xlabel('time (s)');
 ylabel('Dynamin Copies in Cluster');

%Plot experimental curves, scaled to initial density and final relative increase
 plot(edat(:,1)+4,(edat(:,2)-del1)/(edat(1,2)-del1)*acPre(1),'-','Color',[0.5 0.5 0.5], 'LineWidth',3);
 plot(edat(:,1)+4,(edat(:,2)-del2)/(edat(1,2)-del2)*acPost(1),'k-','LineWidth',3);
xlim([0 4])

%Plot model solutions
plot(relPre(:,1),relPre(:,2),'r-','LineWidth',3)
plot(relPost(:,1),relPost(:,2),'m-','LineWidth',3)
