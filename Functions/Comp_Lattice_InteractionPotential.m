function V = Comp_Lattice_InteractionPotential(parm,vars)

% Compute lattice interaction potential: V
rm = 1 + vars.u - vars.u(parm.idx_m);
rp = 1 + vars.u(parm.idx_p) - vars.u;
Vm = Comp_Lattice_LennardJones(parm,rm);
Vp = Comp_Lattice_LennardJones(parm,rp);
V = (Vm+Vp)/2 + parm.V0;

end