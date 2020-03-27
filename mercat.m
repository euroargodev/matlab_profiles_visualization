 function [Axis_factor, Scale_factor] = mercat(lons,lats)%MERCAT  mercator scaling factors for axes 'AspectRatio' parameter.%  This function will calculate scale factors to be used in the axes%  command option 'AspectRatio' to produce mercator projection plots.%  The mercator scaling factor is calculated using Bowditch's formula.%%  The format is as follows:%  [Axis_factor, Scale_factor] = mercat(lons,lats)%%  input:%         lons is a 2 element vector - [minimum lon, maximum lon]%	  lats is a 2 element vector - [minimum lat, maximum lat]%%  output:%         Axis_factor  - scalar%         Scale_factor - scalar%%  usage example:%         axes(... ,'AspectRatio',[Axis_factor,Scale_factor], ...)%%  See the "axes" command discussion in the MATLAB Reference Guide%  for more information on AspectRatio.%% 	Mike Cook - NPS Oceanography Dept. -    v1.0 Oct 93%%   Added error statements                 -    v1.1 Mar 94%   Changed mercator calculation from%   cosine of midpoint latitude to Bowditch%   formula.				   -    v1.2 JUN 94 if nargin ~= 2    error(' Must supply 2 element longitude and latitude vectors ... type help mercat') end if length(lons) ~= 2  |  length(lats) ~= 2    error(' Must supply 2 element longitude and latitude vectors ... type help mercat') end lat_mid = (lats(2) + lats(1)) / 2; L = lat_mid*pi/180;%  M = distance in meters of 1 degree of latitude%  P = distance in meters of 1 degree of longitude%  (source: Bowditch, "The American Practical Navigator")% M = 111132.09 - 566.05*cos(2*L) + 1.2*cos(4*L)- 0.002*cos(6*L); P = 111415.13*cos(L) - 94.55*cos(3*L) + 0.12*cos(5*L); Scale_factor = P/M; Axis_factor = ( (lons(2) - lons(1)) / (lats(2) - lats(1)) ) * Scale_factor;