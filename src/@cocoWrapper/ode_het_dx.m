function J = ode_het_dx(obj, t, z, p, data)
% ODE_HET This function presents nonvectorized implementation of vector field
% Input z - state, p - [omega,epsilon]

assert(~isempty(data),'data to ode_het is empty');

n = obj.system.n;
nt = size(z,2);
% xd = z(n+1:2*n,:);
% om = p(1:end-1,:);
% ep = p(end,:);

% linear part
J1 = repmat(full(obj.system.BinvA),[1,1,nt]);

% nonlinear part
if ~isempty(data) && isfield(data,'fnl')
    fnl = data.fnl; % fnl.coeffs and fnl.ind
    numNonlinearTerms = size(fnl.coeffs,2);
    J2 = zeros(n,2*n,nt);
    Minv = obj.system.M\speye(n);
    for i=1:numNonlinearTerms
        coeff = fnl.coeffs(:,i);
        coeff = Minv*coeff;
        ind = fnl.ind(i,:)';
        % find nonzero exponents
        expind = find(ind);
        xind   = z(expind,:);
        xind(xind==0) = eps; % handle divided by zero
        s = prod(xind.^ind(expind),1);    
        s = (ind(expind)./xind).*s; % #x X nt
        s = s'; s = s(:); s = s';
        s = coeff.*s; s = reshape(full(s),[n,nt,numel(expind)]); s = permute(s,[1,3,2]);
    %     J2(:,ind(expind),:) = J2(:,ind(expind),:)+ s;
        J2(:,expind,:) = J2(:,expind,:)+ s;
    end
else
    J2 =[obj.system.compute_dfnldx(x,xd), obj.system.compute_dfnldxd(x,xd)];
end

J = J1 + [zeros(n,2*n,nt); -J2];

% numerical difference
% f  = @(x,p) obj.ode_het(p(1,:), x, p(2:end,:), data);
% tp = [ t ; p ];
% Jx = coco_ezDFDX('f(x,p)v', f, z, tp);
% dJ = J-Jx;
% norm(dJ(:),'inf')/norm(J(:),'inf')

end
