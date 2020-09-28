folder = fileparts(mfilename('fullpath'));
%% Add GetFullPath folder to MATLAB path
% see: https://www.mathworks.com/matlabcentral/fileexchange/28249-getfullpath
% see: http://blogs.mathworks.com/pick/2011/04/01/be-absolute-about-your-relative-path-with-getfullpath/
addpath(fullfile(folder,'__GetFullPath_20190528'),'-end');

%% Add Units folder to MATLAB path
% see: https://www.mathworks.com/matlabcentral/fileexchange/38977-physical-units-toolbox
addpath(fullfile(folder,'__Sky-s-Physical-Units-for-Matlab-7385215'),'-end');

%% Add export_fig folder to MATLAB path
% see: https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig
addpath(fullfile(folder,'__Altmany-export_fig-9974724'),'-end');

%% Save MATLAB path for successive work sessione
savepath