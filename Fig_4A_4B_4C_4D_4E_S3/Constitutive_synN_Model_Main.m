%%% Main parameter estimation program and manuscript figures

clear all
warning('off')

global divOm_Ha divOs_Ha divOm_La divOs_La yF dat_scal_ConH dat_scal_ConL dat_scal_synH dat_scal_synL
global divOm_Ha_syn divOs_Ha_syn divOm_La_syn divOs_La_syn
global cutoff_time simTime cal lambR  calR cal2 calR2 minS binS maxS

%%%%%%%%%%%%%%%%%% Constitutive
cal2=58;                 %H data converion factor       
calR2=1.8; divRx=calR2;  %R data converion factor      

%%%%%%%%%%%%%%%%%% synN
cal=58;                  %H data converion factor 
calR=1.8;  divRy=calR;   %R data converion factor 

simTime=3;               % Time length (No. of days) of model simulation 
cutoff_time=simTime;     % Same simulation time is used as the cutoff time for the analysis of the data
lambR=0;                 % T cell death rate =0
%%%%%%%

% initial value of parameters (Left: notations as per the model in the manuscript) 
k =2.0783e-08;                % λc -coefficient of target cells lysis   
Kph =0.0072;                  % kp -Transition rate of states in the kinetic proof reading model 
N_hd =7.3789;                 % N  -No. of states formed in the kinetic proof reading model 
A =5.3678e-09;                % ρc -coefficint of T cell proliferation
muC =7.2175;                  % μc -mean of initial T cell distrobution w.r.t. R -Constitutive model
sigC =0.5693;                 % σc -stdv of initial T cell distrobution w.r.t. R -Constitutive model
muS = 6.4;                    % μs of μs.(mean(H)^nH/(mean(H)^nH +KH^nH)) assumed to be the mean of the initial LN T cell distrobution w.r.t. R -synN Model
sigS =1.1;                    % σs -stdv of initial LN T cell distrobution w.r.t. R -synN model
hil =2.4182e+05;              % KH -defletion point of the Hill function
np =3.9704;                   % nH -Hill coefficint

x0=sqrt([k, Kph, N_hd, A, muC, sigC, muS, sigS, hil, np]); % Initial value array
minS=250;binS=500;maxS=8000; % starting bin level, bin sizes, and ending bin level of initial LN T cell distributions w.r.t R    
rng(914052671); % for repeatability
%%%%%%%%%%%%%

%%%%%%%%%%%%% %lysis of target cell data Vs. mean log10(H) from distinct experimentations from Hernadez & Lopez et al. (2021) on day 3:
% Const_High_Affinity
dat_cHD=([3.20	44.79; 4.77	97.75; 5.39	88.13; 5.59	98.46; 6.29	98.10; 7.00	94.45]);
% Const_Low_Affinity
dat_cLD=([3.20	0.34; 4.80	69.12; 5.40	67.67; 5.60	70.37; 6.30	80.10; 7.00	88.18]);
% synN_High_Affinity
dat_sHD=([3.20	0.45; 4.77	30.56; 5.39	29.38; 5.59	69.88; 6.29	92.14; 7.00	79.08]);
% synN_Low_Affinity
dat_sLD=([3.20	0.45; 4.80	3.71; 5.40	0.15; 5.60	12.46; 6.30	46.59; 7.00	58.16]);

%%%%%%%%%%%%% CONSTITUTIVE data:
% Means, stdev, and kurtosys of T cell distributions w.r.t. R in Constitutive experiments on day 3
[meanTR1c,meanTR2c,stdTR1c,stdTR2c,kurtTR1c,kurtTR2c]=meanstd_CONS_Tr(divRx); % statistics from initial Constitutive T cell distribution data files
divOm_Ha=stdTR1c; divOs_Ha=kurtTR1c; divOm_La=stdTR2c;divOs_La=kurtTR2c;

% stdev of % lysis data on day 3 converted from standard error of the mean SEM from graphs: 2.stdev=SEM/sqrt(no.of samples) 
dat_scal_ConH=(sqrt(3).*([12.89, 4.89, 16.97, 5.38, 5.06, 5.06]))./2;                
dat_scal_ConL=(sqrt(3).*([14.11, 4.89, 8.81, 6.69, 4.57, 5.22]))./2;                 
dat_scal_synH=(sqrt(3).*([10.77, 7.13, 5.24, 4.95, 7.86]))./2;                       
dat_scal_synL=(sqrt(3).*([4.66, 5.02, 8.15, 8.15, 5.09]))./2;                        

% Means of target cells U distributions w.r.t. H in Constitutive High Affinity experiments on day3 
MNT1=([meanTR1c,meanTR1c,meanTR1c,meanTR1c,meanTR1c,meanTR1c])';
% Stdev of target cells U distributions w.r.t. H in Constitutive High Affinity experiments on day3
SGT1=([stdTR1c,stdTR1c,stdTR1c,stdTR1c,stdTR1c,stdTR1c])';
SGT1=SGT1.^2; % variance of the same above
% Means of target cells U distributions w.r.t. H in Constitutive Low Affinity experiments on day 3
MNT2=([meanTR2c,meanTR2c,meanTR2c,meanTR2c,meanTR2c,meanTR2c])';
% Stdev of target cells U distributions w.r.t. H in Constitutive High Affinity experiments on day 3
SGT2=([stdTR2c,stdTR2c,stdTR2c,stdTR2c,stdTR2c,stdTR2c])';
SGT2=SGT2.^2; % variance of the same above

% Weighted %lysis of target cells, and mean and variance of T cell distributions w,r,t, R on day 3 -High Affinity, 
my1c=[dat_cHD(1:end,2)./dat_scal_ConH';[MNT1./divOm_Ha;SGT1./divOs_Ha]];
% Weighted %lysis of target cells, and mean and variance of T cell distributions w,r,t, R on day 3 -Low Affinity, 
my2c=[dat_cLD(1:end,2)./dat_scal_ConL';[MNT2./divOm_La;SGT2./divOs_La]];
xc=[my1c,my2c]; yc=([my1c,my2c]);

%%%%%%%%%%%%% synN data:
% Means, stdev, and kurtosys of T cell distributions w.r.t. R in synN experiments on day 3
[meanTR1,meanTR2,meanTR3,meanTR4,meanTR5,meanTR6,stdTR1,stdTR2,stdTR3,stdTR4,stdTR5,stdTR6,...
    kurtTR1,kurtTR2,kurtTR3,kurtTR4,kurtTR5,kurtTR6]=meanstd_SYN_Tr(divRy); % statistics extraction from initial synN T cell distribution data files

divOm_Ha_syn=([stdTR2,stdTR3,stdTR4,stdTR5,stdTR6])'; % stdev of T cell distributions w.r.t. R on day 3 -High Affinity  
divOs_Ha_syn=([kurtTR2,kurtTR3,kurtTR4,kurtTR5,kurtTR6])'; % kurtosys of T cell distributions w.r.t. R on day 3 -High Affinity
divOm_La_syn=([stdTR2,stdTR3,stdTR4,stdTR5,stdTR6])'; % stdev of T cell distributions w.r.t. R on day 3 -Low Affinity
divOs_La_syn=([kurtTR2,kurtTR3,kurtTR4,kurtTR5,kurtTR6])'; % kurtosys of T cell distributions w.r.t. R on day 3 -High Affinity

% Means of target cells U distributions w.r.t. H in synN High & Low Affinity experiments on day3 
MNT=([meanTR2,meanTR3,meanTR4,meanTR5,meanTR6])';
% Stdev of target cells U distributions w.r.t. H in synN High & Low Affinity experiments on day3
SGT=([stdTR2,stdTR3,stdTR4,stdTR5,stdTR6])';
SGT=SGT.^2; % variance of the same above

% Weighted %lysis of target cells, and mean and variance of T cell distributions w,r,t, R on day 3 -High Affinity, 
my1s=[dat_sHD(2:end,2)./dat_scal_synH';[MNT./divOm_Ha_syn;SGT./divOs_Ha_syn]];
% Weighted %lysis of target cells, and mean and variance of T cell distributions w,r,t, R on day 3 -Low Affinity, 
my2s=[dat_sLD(2:end,2)./dat_scal_synL';[MNT./divOm_La_syn;SGT./divOs_La_syn]];

%%%%%%%%%%%%%%%
% Data handling to insert to functions 
xs=[my1s,my2s];
ys=[my1s,my2s];
x=([xc]);
y=([yc]);
yD=([my1c;my2c;my1s;my2s])'; 
xD=([my1c;my2c;my1s;my2s])';

% Non linear fit function
[x00,residual,jacobian,covb]=nlinfit(xD,yD,@Objectivefn_Main_CI_deathTf_estconv_ODE45__CONS_Fin,x0);
% Confidence intervals display
C = nlparci(x00,residual,'covar',covb, 'alpha', 0.05).^2 %CI
% Residual sum of squares
RSS=sum(sum((residual.^2))) % Residual Sum of Squares display
residuals=residual.^2; % Residual squares vector
x00_st=x00; % Storing parameter values 

%%%%%%%% Parameter value display
x00=x00_st;
k=x00(1)^2
Kph=x00(2)^2
N_hd=x00(3)^2
A=x00(4)^2
muC=x00(5)^2
sigC=x00(6)^2
muS=x00(7)^2
sigS=x00(8)^2
hil=x00(9)^2
np=x00(10)^2


%%%%%%%%%%%%%% Preparing const function data for plotting after weighting:  
xD_dat_C=([dat_cHD(1:end,2);dat_cLD(1:end,2)])'./[dat_scal_ConH,dat_scal_ConL]; % Weighted %lysis of target cells on day 3 -Constitutive.
xD_dat_C_M=([MNT1./divOm_Ha;MNT2./divOm_La]); % Weighted mean values of T cell distributions w.r.t. R on day 3 -Constitutive. 
xD_dat_C_S=([SGT1./divOs_Ha;SGT2./divOs_La]); % Weighted variances of T cell distributions w.r.t. R on day 3 -Constitutive.
xD_dat_S=([dat_sHD(2:end,2);dat_sLD(2:end,2)])'./[dat_scal_synH,dat_scal_synL]; % Weighted %lysis of target cells on day 3 -synN.
xD_dat_S_M=([MNT./divOm_Ha_syn;MNT./divOm_La_syn]); % Weighted mean values of T cell distributions w.r.t. R on day 3 -synN. 
xD_dat_S_S=([SGT./divOs_Ha_syn;SGT./divOs_La_syn]); % Weighted variances of T cell distributions w.r.t. R on day 3 -synN.

%%%%%%%%%%%%%%% Preparing const function model estimation as same as abve for plotting after weighting:
xD_mod_C=([yF(1:6),yF(19:24)])';
xD_mod_C_M=([yF(7:12),yF(25:30)])';
xD_mod_C_S=([yF(13:18),yF(31:36)])';
xD_mod_S=([yF(37:41),yF(52:56)])';
xD_mod_S_M=([yF(42:46),yF(57:61)])';
xD_mod_S_S=([yF(47:51),yF(62:66)])';


%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting Cost function without weighting 
% Figure 4B
figure(1002); clf; 
xu=linspace(10^-2,10^8);
pp0=loglog(xu,xu,'-','LineWidth',2); hold on
xxD_dat_C=xD_dat_C'.*([dat_scal_ConH,dat_scal_ConL])';   
xxD_mod_C=xD_mod_C.*([dat_scal_ConH,dat_scal_ConL])'; 
xxD_dat_C_M=[xD_dat_C_M(1:6).*divOm_Ha;xD_dat_C_M(7:12).*divOm_La]; 
xxD_mod_C_M=[xD_mod_C_M(1:6).*divOm_Ha;xD_mod_C_M(7:12).*divOm_La]; 
xxD_dat_C_S=[xD_dat_C_S(1:6).*divOs_Ha;xD_dat_C_S(7:12).*divOs_La]; 
xxD_mod_C_S=[xD_mod_C_S(1:6).*divOs_Ha;xD_mod_C_S(7:12).*divOs_La];
xxD_dat_S=xD_dat_S'.*([dat_scal_synH, dat_scal_synL])';   
xxD_mod_S=xD_mod_S.*([dat_scal_synH, dat_scal_synL])'; 
xxD_dat_S_M=[xD_dat_S_M(1:5).*divOm_Ha_syn;xD_dat_C_M(6:10).*divOm_La_syn]; 
xxD_mod_S_M=[xD_mod_S_M(1:5).*divOm_Ha_syn;xD_mod_C_M(6:10).*divOm_La_syn]; 
xxD_dat_S_S=[xD_dat_S_S(1:5).*divOm_Ha_syn;xD_dat_C_S(6:10).*divOm_La_syn]; 
xxD_mod_S_S=[xD_mod_S_S(1:5).*divOm_Ha_syn;xD_mod_C_S(6:10).*divOm_La_syn];
pp1=loglog(xxD_dat_C,   xxD_mod_C,'o', 'MarkerSize',8, 'linewidth',1); 
pp2=loglog(xxD_dat_C_M, xxD_mod_C_M,'+','MarkerSize',8, 'linewidth',1); 
pp3=loglog(xxD_dat_C_S, xxD_mod_C_S,'^','MarkerSize',8, 'linewidth',1);
pp4=loglog(xxD_dat_S,   xxD_mod_S,'square','MarkerSize',8, 'linewidth',1); 
pp5=loglog(xxD_dat_S_M, xxD_mod_S_M,'*','MarkerSize',8, 'linewidth',1); 
pp6=loglog(xxD_dat_S_S, xxD_mod_S_S,'>','MarkerSize',8, 'linewidth',1); 
legend([pp0 pp1(1) pp2(1) pp3(1) pp4(1) pp5(1) pp6(1)],{'1:1','Const %Lysis','Const T dist means','Const T dist var',...
    'synN %Lysis','synN T dist means','synN T dist var.'}, 'Location', 'northwest')
xlabel('Computations from the data','Fontsize',12)
ylabel('Estimations from the model','Fontsize',12)
box on
title('Fig. 4B');

%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting predictions
% Figure 4C
figure(100); clf; hold on 
Lysis_errors=Plot_synN_HD_ODE45_SYNCONS_CONS_FinSH_PLOT_PREDICTION(k, Kph, N_hd, A, muS, sigS, hil, np);


%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %Lysis fit and other 
% Figure 4A

figure(20); clf; hold on

[fig1,xmn_CH, RmnCH, HmnCH]=Plot_synN_HD_ODE45_SYNCONS_CONS_FinCH(k, Kph, N_hd, A, muC, sigC);  % plots of High Affinity Constitutive
[fig2,xmn_SH, RmnSH, HmnSH]=Plot_synN_HD_ODE45_SYNCONS_CONS_FinSH(k, Kph, N_hd, A, muS, sigS, hil, np);  % plots of High Affinity synN
[fig3,xmn_CL, RmnCL, HmnCL]=Plot_synN_LD_ODE45_SYNCONS_CONS_FinCL(k, Kph, N_hd, A, muC, sigC);  % plots of Low Affinity Constitutive
[fig4,xmn_SL, RmnSL, HmnSL]=Plot_synN_LD_ODE45_SYNCONS_CONS_FinSL(k, Kph, N_hd, A, muS, sigS, hil,np);  % plots of Low Affinity synN


%%%%%%%%%%%%%% Upper CI -%Lysis
figure(20); hold on;

x00=sqrt(C(:,2));

k=x00(1)^2;
Kph=x00(2)^2;
N_hd=x00(3)^2;
A=x00(4)^2;
muC=x00(5)^2;
sigC=x00(6)^2;
muS=x00(7)^2;
sigS=x00(8)^2;
hil=x00(9)^2;
np=x00(10)^2;
%Plot_synN_HD_ODE45_SYNCONS_CONS_FinCH_CI(k, Kph, N_hd, A, muC, sigC);  % CI plots of High Affinity Constitutive
%Plot_synN_HD_ODE45_SYNCONS_CONS_FinSH_CI(k, Kph, N_hd, A, muS, sigS, hil, np);  % CI plots of High Affinity synN
%Plot_synN_LD_ODE45_SYNCONS_CONS_FinCL_CI(k, Kph, N_hd, A, muC, sigC);  % CI plots of Low Affinity Constitutive
%Plot_synN_LD_ODE45_SYNCONS_CONS_FinSL_CI(k, Kph, N_hd, A, muS, sigS, hil,np);  % CI plots of Low Affinity synN


%%%%%%%%%%%%%% Lower CI -%Lysis
figure(20); hold on;

x00=sqrt(C(:,1));

k=x00(1)^2;
Kph=x00(2)^2;
N_hd=x00(3)^2;
A=x00(4)^2;
muC=x00(5)^2;
sigC=x00(6)^2;
muS=x00(7)^2;
sigS=x00(8)^2;
hil=x00(9)^2;
np=x00(10)^2;
%Plot_synN_HD_ODE45_SYNCONS_CONS_FinCH_CI(k, Kph, N_hd, A, muC, sigC);  % CI plots of High Affinity Constitutive
%Plot_synN_HD_ODE45_SYNCONS_CONS_FinSH_CI(k, Kph, N_hd, A, muS, sigS, hil, np);  % CI plots of High Affinity synN
%Plot_synN_LD_ODE45_SYNCONS_CONS_FinCL_CI(k, Kph, N_hd, A, muC, sigC);  % CI plots of Low Affinity Constitutive
%Plot_synN_LD_ODE45_SYNCONS_CONS_FinSL_CI(k, Kph, N_hd, A, muS, sigS, hil,np);  % CI plots of Low Affinity synN

figure(20); hold on
legend([fig1,fig2,fig3,fig4],{'Const High Aff','Const Low Aff','synN High Aff', 'synN Low Aff'}...
    ,'Location','northwest')
legend box off
xlabel('HER2/cell (log10)', 'Fontsize',14)
ylabel('%Lysis', 'Fontsize',14)
box on
title('FIG. 4A')


