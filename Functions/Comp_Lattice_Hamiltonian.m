function H = Comp_Lattice_Hamiltonian(parm,vars)

% Compute lattice Hamiltonian: H
E = Comp_Lattice_EnergyDensity(parm,vars);
H = sum(E); 

end