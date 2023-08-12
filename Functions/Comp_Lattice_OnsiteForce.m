function G = Comp_Lattice_OnsiteForce(parm,vars)

% Compute lattice onsite force: G = -U'
G = -2*pi*parm.U0*sin(2*pi*vars.u);

end