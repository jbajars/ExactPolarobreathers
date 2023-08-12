function U = Comp_Lattice_OnsitePotential(parm,vars)

% Compute lattice onsite potential: U
U = parm.U0*(1-cos(2*pi*vars.u));

end