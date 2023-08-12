function E = Comp_Lattice_EnergyDensity(parm,vars)

% Compute lattice energy density: E
U = Comp_Lattice_OnsitePotential(parm,vars);
V = Comp_Lattice_InteractionPotential(parm,vars);
K = (vars.p.^2)/2;
E = K + V + U; 

end