function [f] = tractor_cost_star_mod(U,z0,zf,Np,Ns,parameters,Optimization_opt)
[F] = tractor_cost_GN_grad_mod(U,z0,zf,Np,Ns,parameters,Optimization_opt);


f=F'*F;