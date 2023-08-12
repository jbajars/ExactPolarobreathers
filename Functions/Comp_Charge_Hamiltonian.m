function Hc = Comp_Charge_Hamiltonian(parm,vars)

% Compute charge Hamiltonian: Hc
Ec = Comp_Charge_EnergyDensity(parm,vars);
Hc = sum(Ec); 

end