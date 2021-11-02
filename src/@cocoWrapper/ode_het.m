function y = ode_het(obj, t, z, p, data)
% ODE_HET This function presents nonvectorized implementation of vector field
% Input z - state, p - [omega,epsilon]

if isempty(data)
    obj.system.Omega        = p(1:end-1,:);
    obj.system.fext.epsilon = p(end,:);
    y = obj.system.odefun(t, z);
else
    n = obj.system.n;
    nt = size(z,2);
    om = p(1,:);
    ep = p(2,:);
    th = p(3,:);

    % linear part
    y1 = obj.system.BinvA*z;

    % nonlinear part
    fnl = data.fnl; % fnl.coeffs and fnl.ind
    numNonlinearTerms = size(fnl.coeffs,2);
    y2 = 0;
    for i=1:numNonlinearTerms
        coeff = repmat(fnl.coeffs(:,i),[1, nt]);
        ind = fnl.ind(i,:)';
        % find nonzero exponents
        expind = find(ind);
        s = prod(z(expind,:).^ind(expind),1);    
        s = repmat(s, [n, 1]);    
        y2 = y2+coeff.*s;
    end
    
    % external forcing
    assert(~isempty(obj.system.fext), 'no external forcing');
    fext_coeffs = repmat(obj.system.fext.coeffs(:,1), [1, nt]);
    fext_harm   = repmat(ep.*cos(obj.system.fext.kappas(1)*om.*t+th), [n, 1]);
    y3 = 2*fext_coeffs.*fext_harm;
    
    y = y1 + [zeros(n,nt); obj.system.M\(-y2+y3)];
end
end