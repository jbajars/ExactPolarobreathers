function dV = Comp_Lattice_LennardJones_deriv(parm,r)

% Compute derivative of the Lennard-Jones potential: dV(r)
R = (1./r).^6;
dV = -12*parm.V0*(R.^2 - R)./r;

end