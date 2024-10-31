function dfnl = compute_dfnldx(obj,x,xd)
% COMPUTE_DFNLDX This function computes the Jacobian of  nonlinear internal 
% force with respect to the displacements x in a second-order
% mechanical system. Currently, we do not treat velocity dependent
% nonlinearity.

assert(obj.order == 2, ' dfnldx can only be computed for second-order systems')

dfnl = sparse(obj.n,obj.n);
switch obj.Options.notation
    case 'tensor'
        if obj.Options.velDepNon
            fnl = zeros(obj.N,obj.N);
            for j = 1:length(obj.fnl)
                fnl = fnl + expand_tensor_derivative(obj.fnl{j},[x(:);xd(:)]);
            end 
            dfnl = fnl(1:obj.n,1:obj.n);            
        else
            for j = 1:length(obj.fnl)
                dfnl = dfnl + expand_tensor_derivative(obj.fnl{j},x);
            end
        end

    case 'fun_handle'
        if obj.Options.velDepNon
            dfnl = obj.dfnldx(x,xd); 
        else
            dfnl = obj.dfnldx(x);
        end
end
