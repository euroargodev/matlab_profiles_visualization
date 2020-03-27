% Argo launcher script

% clear everything
clear all
clc

%--------------------------------------------------------------------------
% Force close of figures
set(0,'ShowHiddenHandles','on')
delete(get(0,'Children'))

h =  findobj('type','figure');
if ~isempty(h)
   delete(h) 
end

h =  findobj('type','uifigure');
if ~isempty(h)
   delete(h) 
end
%--------------------------------------------------------------------------

% Current folder, change it if needed
% root = pwd;
%root = '/data/work_giulio/giulio2';
root = 'C:\Programmi\MATLAB\work\float\medargo\matlab\DMQC\Wong\OW_version_1_0\matlab_database_codes';

%float_data = '/home/ncreati/Desktop/data/work_giulio/giulio2/float_data/';
float_data = 'C:\Programmi\MATLAB\work\float\medargo\data\coriolis_all\';

% data folder
%dataPath = fullfile(root, 'data');
dataPath = fullfile(root, 'ctd_puliti_uniti');

% CTD Path
%CTDPath = fullfile(dataPath, 'CTD');
CTDPath = fullfile(dataPath);

coastPath = fullfile(dataPath);

% addpath(root);
% addpath(CTDPath);
% addpath(coastPath);

appPaths = struct();
appPaths.root = root;
appPaths.CTD = CTDPath;
appPaths.coast = coastPath;
appPaths.data = dataPath;
appPaths.float_data = float_data;

% Start the Argo app
main = ArgoUI(appPaths);


