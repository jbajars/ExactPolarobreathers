function G_W = Comp_ChargeLattice_GW_Force(parm,vars)

% Compute charge-lattice off-diagonal force: F = -dHc_du

% Derivatives of the function J
rm = 1 + vars.u - vars.u(parm.idx_m);
rp = 1 + vars.u(parm.idx_p) - vars.u;
dJm = Comp_Charge_J_tau_deriv(parm,rm);
dJp = Comp_Charge_J_tau_deriv(parm,rp);

% Values of a and b
a = vars.a;
ap = vars.a(parm.idx_p);
am = vars.a(parm.idx_m);
b = vars.b;
bp = vars.b(parm.idx_p);
bm = vars.b(parm.idx_m);

% n-1
G_W = dJm.*(am.*a + bm.*b);

% n+1
G_W = G_W - dJp.*(a.*ap + b.*bp);

end