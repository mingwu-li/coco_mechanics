function data = adj_obj_init_data(fdata,optdof)

data.coll_seg  = fdata.coll_seg;
data.ghan      = @ghan;
data.ghan_dx   = @ghan_dx;
data.ghan_dxdx = @ghan_dxdx;
data.optdof    = optdof;

seg  = fdata.coll_seg;
maps = seg.maps;
int  = seg.int;

NCOL = int.NCOL;
NTST = maps.NTST;
xdim = int.dim;
pdim = maps.pdim;

rows = NCOL*NTST*kron(0:(xdim-1), ones(1,xdim));
opt.gdxdxrows1 = repmat(rows, [1 NCOL*NTST]) + ...
  kron(1:NCOL*NTST, ones(1,xdim^2));
cols = reshape(1:xdim*NCOL*NTST, [xdim NCOL*NTST]);
opt.gdxdxcols1 = repmat(cols, [xdim 1]);

% Derivative of (1/2N)*gxcn with respect to xbp:
step = 1+xdim*(0:NCOL-1); % NCOL
step = repmat(step(:), [1 xdim]) + repmat(0:xdim-1, [NCOL 1]); % xdim*NCOL
step = repmat(step(:), [1 xdim*(NCOL+1)]) + ...
  (xdim*NCOL*NTST+2+pdim)*repmat(0:xdim*(NCOL+1)-1, [xdim*NCOL 1]); % xdim^2*NCOL*(NCOL+1)
step = repmat(step(:), [1 NTST]) + ...
  (xdim*NCOL+xdim*(NCOL+1)*(xdim*NCOL*NTST+2+pdim))*...
  repmat(0:NTST-1, [xdim^2*NCOL*(NCOL+1) 1]); % xdim^2*NCOL*(NCOL+1)*NTST
opt.gdxdxcols2 = step(:);
opt.gdxdxrows2 = ones(xdim^2*NCOL*(NCOL+1)*NTST, 1);

step = 1:NCOL; % NCOL
step = repmat(step(:), [1 xdim]) + NTST*NCOL*repmat(0:xdim-1, [NCOL 1]); % xdim*NCOL
step = repmat(step(:), [1 xdim*(NCOL+1)]) + ...
  xdim*NTST*NCOL*repmat(0:xdim*(NCOL+1)-1, [xdim*NCOL 1]); % xdim^2*NCOL*(NCOL+1)
step = repmat(step(:), [1 NTST]) + (NCOL+xdim^2*NTST*NCOL*(NCOL+1))*...
  repmat(0:NTST-1, [xdim^2*NCOL*(NCOL+1) 1]); % xdim^2*NCOL*(NCOL+1)*NTST
opt.gdxdxidx = step(:);

% Derivative of (T/2N)*gxcn with respect to T:
% opt.gdxdTrows = ones(xdim*NTST*NCOL, 1);
% opt.gdxdTcols = (xdim*NCOL*NTST+2+pdim)*xdim*(NCOL+1)*NTST + ...
%   (1:xdim*NTST*NCOL)';

% Derivative of (1/2N)*w*g' with respect to xbp:
opt.gdTdxcols = NTST*NCOL*xdim+2 + ...
  (xdim*NTST*NCOL+2+pdim)*(0:xdim*(NCOL+1)*NTST-1)';
opt.gdTdxrows = ones(xdim*(NCOL+1)*NTST, 1);

opt.dJrows = 1;
opt.dJcols = (xdim*NTST*NCOL+2+pdim)*(xdim*NTST*(NCOL+1)+2+pdim);

data.coll_opt = opt;

data = coco_func_data(data);

end
