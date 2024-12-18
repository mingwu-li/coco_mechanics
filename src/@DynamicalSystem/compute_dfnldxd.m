function dfnl = compute_dfnldxd(obj,x,xd)
% COMPUTE_DFNLDXD This function computes the Jacobian of  nonlinear internal 
% force f with respect to the velocities xd in a second-order
% mechanical system. Currently, we do not treat velocity dependent
% nonlinearity.

assert(obj.order == 2, ' dfnldxd can only be computed for second-order systems')

dfnl = sparse(obj.n,obj.n);

switch obj.Options.notation
    case 'tensor'
        if obj.Options.velDepNon
            fnl = zeros(obj.N,obj.N);
            for j = 1:length(obj.fnl)
                fnl = fnl + expand_tensor_derivative(obj.fnl{j},[x(:);xd(:)]);
            end 
            dfnl = fnl(1:obj.n,1+obj.n:2*obj.n);            
        end

    case 'fun_handle'
        if obj.Options.velDepNon
            dfnl = obj.dfnldxd(x,xd);
        end
end

