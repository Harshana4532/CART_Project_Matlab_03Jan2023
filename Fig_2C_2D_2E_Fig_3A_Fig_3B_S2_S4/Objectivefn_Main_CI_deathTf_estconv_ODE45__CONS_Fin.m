function [yF]=Objectivefn_Main_CI_deathTf_estconv_ODE45__CONS_Fin(x0,xd)

global yF 

[yF1]=objective_model_all_synN_HD_CI_deathTf_ODE45_SYNCONS_CONS_FinCH(x0);
[yF2]=objective_model_all_synN_LD_CI_deathTf_ODE45_SYNCONS_CONS_FinCL(x0);

yF=max(min(reshape([yF1,yF2],1,[]),100000),0);

