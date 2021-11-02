function J = ghan_dx(x,varargin)

% x1 = x(1,:);
% x2 = x(2,:);
if isempty(varargin)
    optdof = 1:size(x,1);
else
    optdof = varargin{1};
end
[dim,nt] = size(x);
J  = zeros(1,dim,nt);
J(1,optdof,:) = 2*x(optdof,:);

end
