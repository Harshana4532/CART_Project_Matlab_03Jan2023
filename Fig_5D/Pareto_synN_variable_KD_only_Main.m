% Pareto frontier construction between %Lysis of Healthy cells vs 1/%Lysis of Tumor cells, 
% given the synN model compartment of the complex model, 
% with variable koff  
% and fixed KD values 

clear all
global YY_LD YY_HD cutoff_time

%%%%%%%%%%%%%%%%%% Const
cal2=58;                 %H data converion factor       
calR2=1.8; divRx=calR2;  %R data converion factor      
%%%%%%%%%%%%%%%%%% synN
cal=58;                  %H data converion factor 
calR=1.8; divRy=calR;    %R data converion factor 

% Time length of model simulations
simTime=5; % days
% Point in time used in the analysis 
cutoff_time=simTime; 

%%%%%%%
% Defining %Lysis of Healthy cells & 1/%Lysis of Tumor cells
YY_LD=1; YY_HD=1; 

%%%%%% Best estimated Parameters from the complex model
k =2.0787e-08;      % λc in the model
Kph =0.0072;        % kp in the model
N_hd =7.3735;       % N in the model
A =5.3725e-09;      % ρc in the model
muC =7.2431;        % μc in the model
sigC =0.5691;       % σc in the model
muS =6.4049;        % μs in the model 
sigS =1.0999;       % σs in the model
hil =2.4179e+05;    % KH in the model
np =3.9683;         % nH in the model

%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=muS;
sig=sigS;
%%%%%%%%%%%%%%%%%%%%%%%%%%

% For binning of initial T cell Log Normal distributions w.r.t. R 
minS=250;binS=500;maxS=8000;     

global simTime cal calR cal2 calR2 minS binS maxS 

global k Kph N_hd A mu sig KD hil np

% Pareto options in the objective function 
options = optimoptions('gamultiobj','PopulationSize',100,...
          'ParetoFraction',0.7,'PlotFcn',@gaplotpareto);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The choice of KD values in the simulations 
KD=2.39e+0;    %17.6*(0.6*113*0.002)
%KD=2.39e+4; 
%KD=2.39e+5; 
%KD=2.39e+6; 
%KD=2.39e+7; 
%%%%%%%
%Extra simulations: 
%KD=2.39e+2; 
%%%%%%%
%KD=2.39e+3; 
%KD=(2.39+5)*1e+3; 
%%%%%%%
%KD=(2.39+5)*1e+4; 
%KD=(2.39+2.5)*1e+4; 
%KD=(2.39+7)*1e+4; 
%KD=(2.39+9)*1e+4; KD
%KD=(2.39+15)*1e+4;  
%%%%%%%%
%KD=(2.39+5)*1e+5; 
%KD=(2.39+2.5)*1e+5; 
%KD=(2.39+1)*1e+5; 
%KD=(2.39+7)*1e+5; 
%KD=(2.39+9)*1e+5; 
%KD=(2.39+15)*1e+5; 
%%%%%%%%
%KD=(2.39+5)*1e+6; 
%KD=(2.39+2.5)*1e+6; 
%KD=(2.39+1)*1e+6; 
%KD=(2.39+7)*1e+6; 
%KD=(2.39+9)*1e+6; 
%KD=(2.39+15)*1e+6; 
%KD=(2.39+2)*1e+7; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%koff:
koff_LB=1.0E-5;
koff_UB=1.0E-2;
%koff bounds:
lb=[koff_LB]; ub=[koff_UB];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pareto frontier construction function
[solution,ObjectiveValue] = gamultiobj(@mymulti1_synN_variable_KD_only_13dec2022,1,...
                          [],[],[],[],lb,ub,options);


% Plot
figure(1); clf; hold on
sc2=scatter(ObjectiveValue(1:(length(solution)-1),1),ObjectiveValue(1:(length(solution)-1),2),100,(log10(solution(1:(length(solution)-1),1))),'filled');
colorbar
axis([0,100,0.01,0.04]);
colorbar 
sc2.Marker='o';
grid on;
box on
xlabel('%Lysis Healthy')
ylabel('1/%Lysis Tumor')
title("Fig. 4C")

