function G = Comp_ChargeLattice_G_Force(parm,vars)

% Compute total charge-lattice force: F = -dHc_du
G_D = Comp_ChargeLattice_GD_Force(parm,vars);
G_W = Comp_ChargeLattice_GW_Force(parm,vars);

G = G_D + G_W;

end