function J = ode_het_dt(obj, t, z, p, data)
% ODE_HET This function presents nonvectorized implementation of vector field
% Input z - state, p - [omega,epsilon]

assert(~isempty(data),'data to ode_het is empty');

n = obj.system.n;
nt = size(z,2);
om = p(1,:);
ep = p(2,:);
th = p(3,:);

J = zeros(2*n,nt);
% external forcing
assert(~isempty(obj.system.fext), 'no external forcing');
fext_coeffs = obj.system.fext.coeffs(:,1);
Jt = -2*fext_coeffs.*(ep.*obj.system.fext.kappas(1).*om.*sin(obj.system.fext.kappas(1)*om.*t+th));
J(n+1:2*n,:) = obj.system.M\Jt;

% finite difference
% [m, n] = size(z);
% f  = @(x,p) obj.ode_het(x(1,:), p(1:m,:), p(m+1:end,:), data);
% xp = [ z ; p ];
% Jt = reshape(coco_ezDFDX('f(x,p)v', f, t, xp), [m n]);
% dJ = J-Jt;
% norm(dJ(:),'inf')
end