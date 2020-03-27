function File = getfname(Path)
% GETFNAME  - Get the name of all specified files in a directory
%
% Use as: Fname = getfname(Path)
%
% Inputs: Path = Can Include a directory path and filename filter
%                (Default = './*.m')
%
% Output: Fname  = Character array of filenames

% B. Schlining
% 10 Jul 97


if nargin < 1
   Path = '*.m';
end

FDat  = dir(Path);
[r c] = size(FDat);

if r == 0
   File = [];
else
   File = char(FDat.name);
end