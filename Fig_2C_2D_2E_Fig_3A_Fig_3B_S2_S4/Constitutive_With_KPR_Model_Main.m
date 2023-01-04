%%% Constitutive with KPR model parameter estimation program and manuscript figures

clear all
warning('off')

global divOm_Ha divOs_Ha divOm_La divOs_La yF dat_scal_ConH dat_scal_ConL 
global cutoff_time simTime cal calR cal2 calR2 minS binS maxS
global yF 


%%%%%%%%%%%%%%%%%% Const
cal2=58;                 %H data converion factor       
calR2=1.8; divRx=calR2;  %R data converion factor      

%%%%%%%%%%%%%%%%%% synN
cal=58;                 %H data converion factor 
calR=1.8;  divRy=calR;  %R data converion factor 

% Time length of model simulations
simTime=3; %4; % days
cutoff_time=simTime; 


% Parameter values (initial)  
k =2.0783e-08;
Kph =0.0072;
N_hd =7.3789; 7.2; 5.8627;
A =5.3678e-09;
muC =7.2175;
sigC =0.5693;
muS = 6.4;
sigS =1.1;
hil =2.4182e+05;
np =3.9704;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0=sqrt([k, Kph, N_hd, A, muC, sigC]);% Initial value array
minS=250;binS=500;maxS=10000; % bins for initial LN T cell distributions w.r.t R    
rng(914052671); % for repeatability
%%%%%%%%

%%%%%%%%%%%% %Lysis Data:
% Const_HD
dat_cHD=([3.20	44.79; 4.77	97.75; 5.39	88.13; 5.59	98.46; 6.29	98.10; 7.00	94.45]);
% Const_LD
dat_cLD=([3.20	0.34; 4.80	69.12; 5.40	67.67; 5.60	70.37; 6.30	80.10; 7.00	88.18]);

% CONST
[meanTR1c,meanTR2c,stdTR1c,stdTR2c,kurtTR1c,kurtTR2c]=meanstd_CONS_Tr(divRx);
divOm_Ha=stdTR1c; divOs_Ha=kurtTR1c; divOm_La=stdTR2c;divOs_La=kurtTR2c;

%%%%%%%
dat_scal_ConH=(sqrt(3).*([12.89, 4.89, 16.97, 5.38, 5.06, 5.06]))./2;                 
dat_scal_ConL=(sqrt(3).*([14.11, 4.89, 8.81, 6.69, 4.57, 5.22]))./2;                  



MNT1=([meanTR1c,meanTR1c,meanTR1c,meanTR1c,meanTR1c,meanTR1c])';
SGT1=([stdTR1c,stdTR1c,stdTR1c,stdTR1c,stdTR1c,stdTR1c])';
SGT1=SGT1.^2;
MNT2=([meanTR2c,meanTR2c,meanTR2c,meanTR2c,meanTR2c,meanTR2c])';
SGT2=([stdTR2c,stdTR2c,stdTR2c,stdTR2c,stdTR2c,stdTR2c])';
SGT2=SGT2.^2;

my1c=[dat_cHD(1:end,2)./dat_scal_ConH';[MNT1./divOm_Ha;SGT1./divOs_Ha]];
my2c=[dat_cLD(1:end,2)./dat_scal_ConL';[MNT2./divOm_La;SGT2./divOs_La]];

xc=[my1c,my2c];
yc=([my1c,my2c]);


x=([xc]);
y=([yc]);

yD=([my1c;my2c])'; 
xD=([my1c;my2c])';


xD_dat_C=([dat_cHD(1:end,2);dat_cLD(1:end,2)])'./[dat_scal_ConH,dat_scal_ConL];
xD_dat_C_M=([MNT1./divOm_Ha;MNT2./divOm_La]);
xD_dat_C_S=([SGT1./divOs_Ha;SGT2./divOs_La]);

%%%%%%%%%%%%

e_xD_dat_C=xD_dat_C;
e_xD_dat_C_M=xD_dat_C_M;%([yF(7:12),yF(19:24)]).*[dat_scal_ConH,dat_scal_ConL];
%%%%%%%%%%%%


[x00,residual,jacobian,covb]=nlinfit(xD,yD,@Objectivefn_Main_CI_deathTf_estconv_ODE45__CONS_Fin,x0);
C = nlparci(x00,residual,'covar',covb, 'alpha', 0.05).^2; %CI
RSS=sum(sum((residual.^2))); % Residual Sum of Squares
residuals=residual.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x00_st=x00;
k=x00(1)^2
Kph=x00(2)^2
N_hd=x00(3)^2
A=x00(4)^2
muC=x00(5)^2
sigC=x00(6)^2


%%%%%%%%%
xD_mod_C=([yF(1:6),yF(19:24)])';
xD_mod_C_M=([yF(7:12),yF(25:30)])';
xD_mod_C_S=([yF(13:18),yF(31:36)])';


%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 2.D 
figure(1002); clf; 
xu=linspace(10^-2,10^8);
pp0=loglog(xu,xu,'-','LineWidth',2); hold on
xxD_dat_C=xD_dat_C'.*([dat_scal_ConH,dat_scal_ConL])';   
xxD_mod_C=xD_mod_C.*([dat_scal_ConH,dat_scal_ConL])'; 
xxD_dat_C_M=[xD_dat_C_M(1:6).*divOm_Ha;xD_dat_C_M(7:12).*divOm_La]; 
xxD_mod_C_M=[xD_mod_C_M(1:6).*divOm_Ha;xD_mod_C_M(7:12).*divOm_La]; 
xxD_dat_C_S=[xD_dat_C_S(1:6).*divOs_Ha;xD_dat_C_S(7:12).*divOs_La]; 
xxD_mod_C_S=[xD_mod_C_S(1:6).*divOs_Ha;xD_mod_C_S(7:12).*divOs_La];
pp1=loglog(xxD_dat_C,   xxD_mod_C,'o', 'MarkerSize',8, 'linewidth',1); 
pp2=loglog(xxD_dat_C_M, xxD_mod_C_M,'+','MarkerSize',8, 'linewidth',1); 
pp3=loglog(xxD_dat_C_S, xxD_mod_C_S,'^','MarkerSize',8, 'linewidth',1);
legend([pp0 pp1(1) pp2(1) pp3(1)],{'1:1','Const %Lysis','Const T dist means','Const T dist std.'},'Location','northwest')
xlabel('Computations from the data','Fontsize',12)
ylabel('Estimations from the model','Fontsize',12)
box on
title('FIG. 2D')

%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(20); clf;

[fig1, xmn_CH, RmnCH, HmnCH]=Plot_synN_HD_ODE45_SYNCONS_CONS_FinCH2(k, Kph, N_hd, A, muC, sigC);
[fig2, xmn_CL, RmnCL, HmnCL]=Plot_synN_LD_ODE45_SYNCONS_CONS_FinCL2(k, Kph, N_hd, A, muC, sigC);

figure(20); hold on
legend([fig1,fig2],{'Const High Aff','Const Low Aff'}...
    ,'Location','northwest')
legend box off
xlabel('HER2/cell (log10)', 'Fontsize',14)
ylabel('%Lysis', 'Fontsize',14)
box on
title('FIG. 2C')



%%%%%%%%%%%%%%% Figure S.4
figure(700); clf; hold on
%High Aff Const High -Normal cells
Rmn1=RmnCH(2);
Hmn1=HmnCH(2);
koff_hd=9E-5; % High Aff
coun=1;
for u=0:0.1:8
    KD=17.6.*(10^u).*(0.6*113*0.002);
    bet=1+Kph/koff_hd; KDp=(KD+(Kph/(koff_hd/KD)))/bet; 
    xmnCN_Norm(coun)=0.5.*(Rmn1+Hmn1+KDp).*(1-sqrt(1-(4.*Rmn1.*Hmn1)./((Rmn1+Hmn1+KDp).^2))).*(Kph./(Kph+koff_hd)).^N_hd;
KDCN_Norm(coun)=KD;
coun=coun+1;
end
pl1=plot(log10(KDCN_Norm),xmnCN_Norm,'o');
plot(log10(KDCN_Norm),xmnCN_Norm,'-')
%%%%%%%%%%%%%%%%
%High Aff Const High -Tumor cells
Rmn1=RmnCH(5);
Hmn1=HmnCH(5);
koff_hd=9E-5; % High Aff
coun=1;
for u=0:0.1:8
    KD=17.6.*(10^u).*(0.6*113*0.002);
    bet=1+Kph/koff_hd; KDp=(KD+(Kph/(koff_hd/KD)))/bet; 
    xmnCN_Tum(coun)=0.5.*(Rmn1+Hmn1+KDp).*(1-sqrt(1-(4.*Rmn1.*Hmn1)./((Rmn1+Hmn1+KDp).^2))).*(Kph./(Kph+koff_hd)).^N_hd;
KDCN_Tum(coun)=KD;
coun=coun+1;
end
pl2=plot(log10(KDCN_Tum),xmnCN_Tum,'^');
plot(log10(KDCN_Tum),xmnCN_Tum,'-')
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
%High Aff Const High -Normal cells
Rmn1=RmnCH(2);
Hmn1=HmnCH(2);
koff_hd=9E-5; % High Aff
coun=1;
for u=0:0.1:8
    KD=17.6.*(10^u).*(0.6*113*0.002);
    xmnC_Norm(coun)=0.5.*(Rmn1+Hmn1+KD).*(1-sqrt(1-(4.*Rmn1.*Hmn1)./((Rmn1+Hmn1+KD).^2)));
KDC_Norm(coun)=KD;
coun=coun+1;
end
pl3=plot(log10(KDC_Norm),xmnC_Norm,'square');
plot(log10(KDC_Norm),xmnC_Norm,'-')
%%%%%%%%%%%%%%%%
%High Aff Const High -Tumor cells
Rmn1=RmnCH(5);
Hmn1=HmnCH(5);
koff_hd=9E-5; % High Aff
coun=1;
for u=0:0.1:8
    KD=17.6.*(10^u).*(0.6*113*0.002);
    xmnC_Tum(coun)=0.5.*(Rmn1+Hmn1+KD).*(1-sqrt(1-(4.*Rmn1.*Hmn1)./((Rmn1+Hmn1+KD).^2)));
KDC_Tum(coun)=KD;
coun=coun+1;
end
pl4=plot(log10(KDC_Tum),xmnC_Tum,'>-');
plot(log10(KDC_Tum),xmnC_Tum,'-')
%%%%%%%%%%%%%%%%
xlabel('K_D (log10)', 'Fontsize',12)
ylabel('Number of Complexes at Mean R and Mean HER per cell', 'Fontsize',12)
box on
title('Modeled [RH]/cell, Constitutive High Affinity')
legend([pl1 pl2 pl3 pl4],{'C_N-Healthy','C_N-Tumor',... 
    'C_0-Healthy','C_0-Tumor'})
legend box off
title('Fig. S4') 

%%%%%%%%%%

x00=sqrt(C(:,2));
k=x00(1)^2;
Kph=x00(2)^2;
N_hd=x00(3)^2;
A=x00(4)^2;
muC=x00(5)^2;
sigC=x00(6)^2;


%%%%%%%

x00=sqrt(C(:,1));
k=x00(1)^2;
Kph=x00(2)^2;
N_hd=x00(3)^2;
A=x00(4)^2;
muC=x00(5)^2;
sigC=x00(6)^2;


% Figure 2.E
%%%%%%%%%% Complexes
x00=x00_st;
k=x00(1)^2;
Kph=x00(2)^2;
N_hd=x00(3)^2;
A=x00(4)^2;
muC=x00(5)^2;
sigC=x00(6)^2;

KD=17.6*(0.6*113*0.002);
koff_hd=9E-5;
Rmn1=2842.6;
clear xmnx Hmn1
%Rmn1(p)=sum(R.*Tr(1,:)')./sum(Tr(1,:));
vale=3;
for i=1:40;
Hmn1(i)=10^vale;
    bet=1+Kph/koff_hd; KDp=(KD+(Kph/(koff_hd/KD)))/bet; 
    xmnx(i)=0.5.*(Rmn1+Hmn1(i)+KDp).*(1-sqrt(1-(4.*Rmn1.*Hmn1(i))./((Rmn1+Hmn1(i)+KDp).^2))).*(Kph./(Kph+koff_hd)).^N_hd;
vale=vale+0.1;
end

figure(85); clf
plot(log10(Hmn1),xmnx,'o-'); hold on
    %%%%
KD=210*(0.6*113*0.002);
koff_hd=6.8E-4;
Rmn1=2842.6;
clear xmnx Hmn1
vale=3;
for i=1:40;
Hmn1(i)=10^vale;
    bet=1+Kph/koff_hd; KDp=(KD+(Kph/(koff_hd/KD)))/bet; 
    xmnx(i)=0.5.*(Rmn1+Hmn1(i)+KDp).*(1-sqrt(1-(4.*Rmn1.*Hmn1(i))./((Rmn1+Hmn1(i)+KDp).^2))).*(Kph./(Kph+koff_hd)).^N_hd;
vale=vale+0.1;
end
hold on; figure(85)
plot(log10(Hmn1),xmnx,'o-'); hold on
title('FIG. 2E')
xlabel('HER2 antigen density (molec/cell)')
ylabel('CAR-HER2')
legend('High Affinity','Low Affinity')
axis([-inf,inf,0,3000])

    %%%%


