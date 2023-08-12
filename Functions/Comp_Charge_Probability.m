function C2_density = Comp_Charge_Probability(parm,vars)

% Compute charge density
C2_density = (vars.a.^2 + vars.b.^2)/2/parm.tau;

end