function [data, y] = objhan(prob, data, u) %#ok<INUSL>

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T    = u(maps.T_idx);
x    = u(maps.xbp_idx);
p    = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.ghan(xcn,data.optdof);
gcn = mesh.gka.*gcn;   % mesh here is related adaptive meshing

y = (0.5/maps.NTST)*mesh.gwt*gcn';

end
