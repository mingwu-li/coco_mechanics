function install
% run this file first to install all external packages  and
% also the main software appropriately

maindir = fileparts(mfilename('fullpath'));
addpath(fullfile(maindir, 'src'));
addpath(fullfile(maindir, 'ext', 'tensor_toolbox'));
addpath(fullfile(maindir, 'ext', 'misc'));

end


