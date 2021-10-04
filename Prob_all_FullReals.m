%% Stress on Faults input
load('tableTrim.mat')
% The other Mohr circle codes are in Focal - This one is for simple 2D
% circles of each reservoir%
% PH=13117; %Basement
PH=8200; %Hogback
% PH=7287; %Dakota/Missippian Carbonates/Mesaverde
subn=0;
Reals=2*10^4; %realizations


untrunc = makedist('Normal',0,10);
trunc = truncate(untrunc,0,45);



for i=[1]%,3,5]

    
    
    subn=subn+1;
S1=tableTrim(find(tableTrim(:,1)==PH,1),2);
S2=tableTrim(find(tableTrim(:,1)==PH,1),5)*PH;
S3=tableTrim(find(tableTrim(:,1)==PH,1),4);
Pp=tableTrim(find(tableTrim(:,1)==PH,1),3);
mu=0.6;
a=305;
b=-90;
c=0;
% faults_nodal=[N1_strike N1_dip;
%               N2_strike N2_dip];
%           
% rake_nodal=[N1_rake;
%             N2_rake];
% rake_nodal=rake_nodal(IndexOfNodal);
%           
% faults_nodal=faults_nodal(:,IndexOfNodal);
% 

% Stress on faults code

%converting fault and stress orientation info from degrees to radians 
% str=deg2rad([33,33,33,33,33,33,33,33,33,33,213,213,213,213,213,213,213,213,213]);
%str=deg2rad([100,100,100,100,100,100,100,100,100,100,280,280,280,280,280,280,280,280,280]);
% str=deg2rad([90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 270, 270, 270, 270, 270, 270, 270, 270, 270]);
%dip=deg2rad([0,10,20,30,40,50,60,70,80,90,80,70,60,50,40,30,20,10,0]); % :=1 !!!

% a=deg2rad(a); % deg2rad done at sampling (below)
b=deg2rad(b);
c=deg2rad(c);


for j=1:Reals
% 
% str=deg2rad(normrnd(90,10)); % For Basement and Dakota (35 base, 125 base2, 90 dakota)
% dip=deg2rad(unifrnd(0,90));
% 
% str=deg2rad(normrnd(33,10)); % For Hogback
% dip=deg2rad(normrnd(50,10));
% if dip<0
%     dip=0;
% end

str = deg2rad(normrnd(45,10)); % For Hogback Vertical
dip = deg2rad(random(trunc));




S1_real=normrnd(S1,S1*.04);
S2_real=normrnd(S2,S2*.10);
S3_real=normrnd(S3,S3*.10);
Pp_real=normrnd(Pp,Pp*.04);
mu_real=normrnd(mu,mu*.15);
a_real=deg2rad(normrnd(a,10));

if S2_real < S3_real
    PlaceHold=S2_real;
    S2_real=S3_real;
    S3_real=PlaceHold;
end

% str=deg2rad(33); % For Hogback
% dip=deg2rad(70);
% 
% S1_real=S1;
% S2_real=S2;
% S3_real=S3;
% Pp_real=Pp;
% mu_real=mu;
% a_real=deg2rad(a);

%defines principal stress tensor
S=[S1_real 0 0
    0 S2_real 0
    0 0 S3_real];
%transformation from principal stress to geographic coordinates
R1=[cos(a_real)*cos(b) sin(a_real)*cos(b) -sin(b);
    cos(a_real)*sin(b)*sin(c)-sin(a_real)*cos(c) sin(a_real)*sin(b)*sin(c)+cos(a_real)*cos(c) cos(b)*sin(c);
    cos(a_real)*sin(b)*cos(c)+sin(a_real)*sin(c) sin(a_real)*sin(b)*cos(c)-cos(a_real)*sin(c) cos(b)*cos(c)];

Sg=R1'*S*R1;


%the following part of the code is dependent on fault orientation
%for loop makes calculations for each fault one at a time
% fault_no=[1:1:size(faults_nodal,2)]';
for k=1%:19;%length(fault_no)
    % transformation from geographic to fault coordinate system
    R2=[cos(str(k)) sin(str(k)) 0;
         sin(str(k))*cos(dip(k)) -cos(str(k))*cos(dip(k)) -sin(dip(k));
         -sin(str(k))*sin(dip(k)) cos(str(k))*sin(dip(k)) -cos(dip(k))];

    Sf=R2*Sg*R2';
    %Solves for normal stress resolved on the fault surface
    Sn(k)=Sf(3,3);

    %Solve for rake of the slip vector
    if (Sf(3,2)>0)
        rake(k)=(atan(Sf(3,2)/(Sf(3,1))));
    elseif (Sf(3,2)<0)&(Sf(3,1)>0)
        rake(k)=pi-atan(Sf(3,2)/(-Sf(3,1)));
    else
        rake(k)=atan((-Sf(3,2))/(-Sf(3,1)))-pi;
    end
    %Solve for shear stress resolved on the fault plane
    R3=[cos(rake(k)) sin(rake(k)) 0;
        -sin(rake(k)) cos(rake(k)) 0;
        0 0 1];

    Sr=R3*Sf*R3';

    tau(k,1)=Sr(3,1);
    %solve for Coulomb failure function (CFF) (proxy for likelihood of slip) of fault
    Sn_eff(k,1)=Sn(k)-Pp_real;
    CFF(k,1)=abs(tau(k)) - mu_real*(Sn_eff(k));
    rake(k)=rad2deg(rake(k));
end

CFF_all(j)=CFF(k,1);
tau_all(j)=tau(k,1);
Sn_eff_all(j)=Sn_eff(k,1);

end
end
% %output data
% %changing signs to make output easier to understand
% % rake=(sign(tau).*rake');
% % tau=abs(tau);
% 
% %FaultNo_Strike_Dip_ShearStress_EffNormalStress_Rake_CFF(:,i)=[fault_no,faults_nodal',tau, Sn_eff, rake, CFF]
% 
% % FailureFunc(:,i)=CFF(:);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % Visualization
% 
% %This part of the code is for visualizing the above results on 
% %a 3D Mohr Diagram
% S1_eff=S1_real-Pp_real;
% S2_eff=S2_real-Pp_real;
% S3_eff=S3_real-Pp_real;
% 
% 
% %%%%Building the 3D diagram
% format compact                    
% format long e                    
% eta=linspace(0,pi,50);
% hold on
% % figure(1)
% % subplot(4,2,2*(subn)-1)
% 
% %Plotting S3-S2 cicle
% x_32 = ((S2_eff-S3_eff)/2)*cos(eta)+((S2_eff-S3_eff)/2+S3_eff);                   % generate x-coordinate
% y_32 = ((S2_eff-S3_eff)/2)*sin(eta);                   
% plot(x_32,y_32,'k','Linewidth',.01);                        
% axis('equal');                    
% % set(gca, 'Xlim', [-1, S1_eff+2], 'YLim', [0, (S1(i)-S3(i))/2*1.1])
% set(gca, 'Xlim', [0, 11000], 'YLim', [0, 4000])
% xlabel('Normal Stress, S_{ij} (psi)')
% 
% ylabel('Shear Stress (psi)')
% % title(strcat('\Delta P = ',num2str(16-i+1),'MPa'))
% hold on
% 
% %Plotting S2-S1 cicle
% x_21 = ((S1_eff-S2_eff)/2)*cos(eta)+((S1_eff-S2_eff)/2+S2_eff);                   % generate x-coordinate
% y_21 = ((S1_eff-S2_eff)/2)*sin(eta);                   
% plot(x_21,y_21,'k','Linewidth',.01);                        
% 
% %Plotting S3-S1 cicle
% x_31 = ((S1_eff-S3_eff)/2)*cos(eta)+((S1_eff-S3_eff)/2+S3_eff);                   % generate x-coordinate
% y_31 = ((S1_eff-S3_eff)/2)*sin(eta);                   
% plot(x_31,y_31,'k','Linewidth',.01);                        
% 
% %plot failure line
% y_failure=[0 1000000];%(S1-S3)/2*2];
% x_failure=y_failure/mu_real;
% plot(x_failure, y_failure, 'r-','Linewidth',.01)
% 
% 
% end
% 
% for j=1:Reals
% hold on
% plot(Sn_eff_all(j), abs(tau_all(j)),'.', 'MarkerSize',20, 'Color',[1-(abs(CFF_all(j))-abs(max(CFF_all)))/max(abs(CFF_all)) 0 0]);
% end
% 
% %%%%Adding results from previous calculations
% % scatter(faults(:,1), faults(:,2), 50, CFF, 'filled')
% % title('Fault Orientations (labeled by input order & colorcoded by CFF)')
% % hold on
% % for i=1:length(fault_no)
% % str=num2str(i);
% % text(faults(i,1)+3, faults(i,2)+3,str,'fontsize',12)
% % end
% % set(gca,'Xlim',[0 360], 'Ylim', [0 90])
% % xlabel('Strike')
% % ylabel('Dip')
% % colorbar
% %set(gca, 'Xlim', [0 180], 'Ylim', [-90 90])
% 
% % if subn==1
% %     xlim([27 52])
% % else
% % xlim([27 52])
% % end
% % 
% % title(strcat('\Delta P = ',num2str(15-i+1),'MPa'))
% % for ii = 1:size(FailureFunc,1)
% % if FailureFunc(ii,i)>0
% % scatter(Sn_eff(ii), abs(tau(ii)), 50, 'r', 'filled')
% % else
% % scatter(Sn_eff(ii), abs(tau(ii)), 50, 'k', 'filled')
% % end
% % format short
% % end
% 
% 
% 
% 
% S1_eff=S1-Pp;
% S2_eff=S2-Pp;
% S3_eff=S3-Pp;
% 
% 
% %%%%Building the 3D diagram
% format compact                    
% format long e                    
% eta=linspace(0,pi,50);
% hold on
% % figure(1)
% % subplot(4,2,2*(subn)-1)
% 
% %Plotting S3-S2 cicle
% x_32 = ((S2_eff-S3_eff)/2)*cos(eta)+((S2_eff-S3_eff)/2+S3_eff);                   % generate x-coordinate
% y_32 = ((S2_eff-S3_eff)/2)*sin(eta);                   
% plot(x_32,y_32,'k');                        
% axis('equal');                    
% % set(gca, 'Xlim', [-1, S1_eff+2], 'YLim', [0, (S1(i)-S3(i))/2*1.1])
% set(gca, 'Xlim', [0, 10500], 'YLim', [0, 3500])
% xlabel('Effective Normal Stress (psi)')
% 
% ylabel('Shear Stress (psi)')
% % title(strcat('\Delta P = ',num2str(16-i+1),'MPa'))
% hold on
% 
% %Plotting S2-S1 cicle
% x_21 = ((S1_eff-S2_eff)/2)*cos(eta)+((S1_eff-S2_eff)/2+S2_eff);                   % generate x-coordinate
% y_21 = ((S1_eff-S2_eff)/2)*sin(eta);                   
% plot(x_21,y_21,'k');                        
% 
% %Plotting S3-S1 cicle
% x_31 = ((S1_eff-S3_eff)/2)*cos(eta)+((S1_eff-S3_eff)/2+S3_eff);                   % generate x-coordinate
% y_31 = ((S1_eff-S3_eff)/2)*sin(eta);                   
% plot(x_31,y_31,'k');                        
% 
% %plot failure line
% y_failure=[0 1000000];%(S1-S3)/2*2];
% x_failure=y_failure/mu_real;
% plot(x_failure, y_failure, 'r-')
% 
% 
% 
% 
% 
% end
