function F = Comp_Lattice_InteractionForce(parm,vars)

% Compute lattice interaction force: F = -V'
rm = 1 + vars.u - vars.u(parm.idx_m);
rp = 1 + vars.u(parm.idx_p) - vars.u;
dVm = Comp_Lattice_LennardJones_deriv(parm,rm);
dVp = Comp_Lattice_LennardJones_deriv(parm,rp);
F = -(dVm - dVp);

end