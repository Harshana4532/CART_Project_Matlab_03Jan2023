function [f]= mymulti1_syn_variable_nH_KH_fixed_KD_13dec2022(x)

global k Kph N_hd A mu sig YY_LD YY_HD

% Variables used in optimization function 
hil=x(1); % K_H Hill function parameter 
np=x(2); % n_H Hill function parameter  

% Model that generates the %Lysis of Healthy (YY_LD) vs % 1/%Lysis of tumor cells (YY_HD)
[YY_LD,YY_HD]=pareto_model_syn_variable_nH_KH_fixed_KD_13dec2022(k, Kph, N_hd, A, mu, sig, np, hil );

% The two objective functions 
f(1) =YY_LD; 
f(2) =YY_HD; 

end

