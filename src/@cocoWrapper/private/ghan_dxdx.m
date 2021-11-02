function J = ghan_dxdx(x,varargin)

[nx,nt] = size(x);
if isempty(varargin)
    optdof = 1:size(x,1);
else
    optdof = varargin{1};
end
J  = zeros(1, nx, nx, nt);

for k=1:numel(optdof)
    J(1,optdof(k),optdof(k),:) = 2;
end

end