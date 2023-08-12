function G_B = Comp_ChargeLattice_GB_Force(parm,vars)

% Compute charge-lattice off-diagonal B force: F = -dHc_du

% Derivatives of the function J
rm = 1 + vars.u - vars.u(parm.idx_m);
rp = 1 + vars.u(parm.idx_p) - vars.u;
dJm = Comp_Charge_J_tau_deriv(parm,rm);
dJp = Comp_Charge_J_tau_deriv(parm,rp);

% Values of a 
a = vars.a;
ap = vars.a(parm.idx_p);
am = vars.a(parm.idx_m);

% n-1
G_B= dJm.*(am.*a);

% n+1
G_B = G_B - dJp.*(a.*ap);

end