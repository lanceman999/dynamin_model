
%time in s
%conc in uM
%V in um3
%Atotal and Acluster in um2
function[relTot]=plot_dynamin_fromODE(time, conc, V, Atotal, Acluster, n)

copies=conc*602*V;%conc is in units of uM. convert to copies.
c=copies(:,4)/Atotal;
ac=copies(:,2)/Atotal;
acd=copies(:,5)/Acluster;
ad=copies(:,6)/Acluster;

f=figure(n);
 ax=axes('Parent',f,'FontSize',20,'LineWidth',1);
 hold(ax);
 xlabel('time (s)');
 ylabel('copies/\mum^2 (on surface)');
%recruiter and dynamin both dilute phase
plot(time, c(:,1),'b--');
plot(time, ac(:,1),'c--');

%dynamin in cluster 
plot(time, acd(:,1),'r--');%from 2D
plot(time, ad(:,1),'m--');%from 3D

legend('recruiter in dilute phase','dynamin in dilute phase','dynamin in cluster (from 2D)','dynamin in cluster (from 3D)')

f2=figure(n+1);
 ax2=axes('Parent',f2,'FontSize',20,'LineWidth',1);
hold(ax2);


xlabel('time (s)');
 ylabel('Dynamin Copies');
 cmap=colormap('jet');
%recruiter and dynamin both dilute phase
%plot(time, c*Atotal,'b--');
plot(time, ac*Atotal,'b-','LineWidth',2);
%dynamin in cluster 
plot(time, acd*Acluster,'c-','LineWidth',2);%from 2D
plot(time, ad*Acluster,'-','LineWidth',2,'Color',cmap(130,:));%from 3D
plot(time, ad*Acluster+acd*Acluster+ac*Acluster,'k-','LineWidth',3);%from 3D
legend('[D]_{mem}','[D]_{clus2D}','[D]_{clus3D}','[D]_{clus} total')

f3=figure(n+2);
 ax3=axes('Parent',f3,'FontSize',20,'LineWidth',1);
hold(ax3);


xlabel('time (s)');
 ylabel('Dynamin Copies in Cluster');
%recruiter and dynamin both dilute phase
display('Copies initially incluster')
ac(1)*Acluster
display('Copies in cluster at end');
acd(end)*Acluster+ad(end)*Acluster+ac(end)*Acluster

display('Copies initially outside of cluster')
ac(1)*Atotal
display('Copies outside of cluster at end')
ac(end)*Atotal


%Intensity is relative to the value at time zero. 
%plot(time, (acd+ad+ac)./ac(1),'k-','LineWidth',1);
%Here is just the change in copies in the cluster
plot(time, (acd+ad+ac)*Acluster,'k-','LineWidth',1)

relTot=[time', (acd+ad+ac)*Acluster];
