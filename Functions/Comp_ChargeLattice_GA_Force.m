function G_A = Comp_ChargeLattice_GA_Force(parm,vars)

% Compute charge-lattice off-diagonal A force: F = -dHc_du

% Derivatives of the function J
rm = 1 + vars.u - vars.u(parm.idx_m);
rp = 1 + vars.u(parm.idx_p) - vars.u;
dJm = Comp_Charge_J_tau_deriv(parm,rm);
dJp = Comp_Charge_J_tau_deriv(parm,rp);

% Values of b
b = vars.b;
bp = vars.b(parm.idx_p);
bm = vars.b(parm.idx_m);

% n-1
G_A = dJm.*(bm.*b);

% n+1
G_A = G_A - dJp.*(b.*bp);

end