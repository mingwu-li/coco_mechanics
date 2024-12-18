function fnl = compute_fnl(obj,x,xd)
% COMPUTE_FNL We compute the nonlinear internal force in a second-order
% mechanical system. Currently, we do not treat velocity dependent
% nonlinearity.

assert(obj.order == 2, ' fnl can only be computed for second-order systems')

switch obj.Options.notation
    case 'tensor'
        if obj.Options.velDepNon
            fnl = zeros(obj.N,1);
            for j = 1:length(obj.fnl)
                fnl = fnl + expand_tensor(obj.fnl{j},[x(:);xd(:)]);
            end 
            fnl = fnl(1:obj.n);
        else
            fnl = zeros(obj.n,1);
            for j = 1:length(obj.fnl)
                fnl = fnl + expand_tensor(obj.fnl{j},x);
            end
        end
    case 'fun_handle'
        if obj.Options.velDepNon
            fnl = obj.fnl(x,xd);
        else
            fnl = obj.fnl(x);
        end
    otherwise
        error('notation should be tensor or fun_handle');
end