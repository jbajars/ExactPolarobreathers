function [H,Hlat,Hc] = Comp_Hamiltonian(parm,vars)

% Compute total Hamiltonian

% Lattice Hamiltonian: Hlat
Hlat = Comp_Lattice_Hamiltonian(parm,vars);

% Charge Hamiltonian: Hc
Hc = Comp_Charge_Hamiltonian(parm,vars);

% Total Hamiltonian: H
H = Hlat + Hc;

end