function [data, y] = adj_objhan(prob, data, u) %#ok<INUSL>

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T = u(maps.T_idx);
x = u(maps.xbp_idx);
p = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);

gdxcn = pr.ghan_dx(xcn,data.optdof);
gdxcn = mesh.gdxka.*gdxcn;

J_xbp = (0.5/maps.NTST)*gdxcn(:)';
J_T0  = 0;
J_T   = 0;
J_p   = zeros(1,numel(p));

y = [J_xbp J_T0 J_T J_p];

end
