function [data, y] = OmegaT(prob, data, u) %#ok<INUSL>
y = u(1)*u(2)-2*pi*data.periodsRatio; % omega*T = 2*pi
end