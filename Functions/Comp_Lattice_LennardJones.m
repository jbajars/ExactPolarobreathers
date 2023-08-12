function V = Comp_Lattice_LennardJones(parm,r)

% Compute Lennard-Jones potential: V(r)
r = (1./r).^6;
V = parm.V0*(r.^2 - 2*r);

end