function Ec = Comp_Charge_EnergyDensity(parm,vars)

% Compute charge energy density: Ec
Ec = parm.E0*Comp_Charge_Probability(parm,vars);

% Compute function J values
um = vars.u - vars.u(parm.idx_m);
up = vars.u(parm.idx_p) - vars.u;
rm = 1 + um;
rp = 1 + up;
Jm = Comp_Charge_J_tau(parm,rm);
Jp = Comp_Charge_J_tau(parm,rp);

% Compute a and b values
a = vars.a;
b = vars.b;
am = vars.a(parm.idx_m);
ap = vars.a(parm.idx_p);
bm = vars.b(parm.idx_m);
bp = vars.b(parm.idx_p);

Ec = Ec - (Jp.*(a.*ap + b.*bp) + Jm.*(a.*am + b.*bm))/2;

end