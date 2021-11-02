function y = ghan(x,varargin)

% x1 = x(1,:);
% x2 = x(2,:);
% 
% y  = x1./(1+x2.^2);
if isempty(varargin)
    optdof = 1:size(x,1);
else
    optdof = varargin{1};
end
x = x(optdof,:);

y = sum(x.^2,1);

end
