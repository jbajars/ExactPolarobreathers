function C2 = Comp_Charge_Total(parm,vars)

% Compute total charge
C2_density = Comp_Charge_Probability(parm,vars);
C2 = sum(C2_density);

end