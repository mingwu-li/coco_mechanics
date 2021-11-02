function J = ode_het_dp(obj, t, z, p, data)
% ODE_HET This function presents nonvectorized implementation of vector field
% Input z - state, p - [omega,epsilon]

assert(~isempty(data),'data to ode_het is empty');

n = obj.system.n;
nt = size(z,2);
om = p(1,:);
ep = p(2,:);
th = p(3,:);

J = zeros(2*n,3,nt);
% external forcing
assert(~isempty(obj.system.fext), 'no external forcing');
fext_coeffs = obj.system.fext.coeffs(:,1);
Jom = -2*fext_coeffs.*(ep.*obj.system.fext.kappas(1).*t.*sin(obj.system.fext.kappas(1)*om.*t+th));
Jep = 2*fext_coeffs.*cos(obj.system.fext.kappas(1)*om.*t+th);
Jth = -2*fext_coeffs.*(ep.*sin(obj.system.fext.kappas(1)*om.*t+th));
J(n+1:2*n,1,:) = reshape(Jom,[n,1,nt]);
J(n+1:2*n,2,:) = reshape(Jep,[n,1,nt]);
J(n+1:2*n,3,:) = reshape(Jth,[n,1,nt]);

% finite difference
% f  = @(x,p) obj.ode_het(x(1,:), x(2:end,:), p, data);
% tx = [ t ; z ];
% Jp = coco_ezDFDP('f(x,p)v', f, tx, p);
% dJ = J-Jp;
% norm(dJ(:),'inf')

end