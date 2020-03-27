% =================================================================================
%
% argo_file_read          reads some of the variables contained in an
%                         ARGO format (profile) file (netcdf). See
%                         reference:
%                         "Argo User's manual" available on
%                         http://www.coriolis.eu.org/coriolis/cdc/argo_rfc.htm
%
% AUTHOR:                Emmanuelle Autret
%
% HISTORIC:              23/10/2003 - EA - version 1.00 - 
%                        with matlab version 6.13
%			 30/08/2004 -TL - version 1.01 -
%			 update for argo format version 2.01

% REQUIRES:              netcdf matlab toolbox
%
% EXAMPLE:               file_name='6900162_500621518.nc';
%
%                        [latitude,longitude,juld,pressure,temperature,salinity,platform_number,cycle_number,data_type,format_version,project_name,pi_name,direction,data_centre,data_state_indicator,data_mode,dc_reference,inst_reference,wmo_inst_type]=argo_file_read(file_name);
%
%
% SEE ALSO:              plot_latest_TS_profile.m, plot_contour.m.      
% 
% ================================================================================
% MODIFIED by G. Notarstefano 13/03/2008 ==> "position_qc, pres_qc, psal_qc, temp_qc"
% added in output
% =========================================================================


function argoStruct = argo(file_name, app);
%  function argoStruct = argo();

%  file_name = '/data/work_giulio/giulio2/float_data/6901655/profiles/D6901655_008.nc'

%ncstartup

% ncf = netcdf.open(file_name,'NOWRITE'); 
% ncf = netcdf(file_name,'read');
% variables = var(ncf);
%Reads genaral informations
%============================

% MULTIPROFILI
% nel caso di multiprofili, considero il PRIMARY SAMPLING che dovrebbe
% sempre essere nella prima colonna

% *************************************************************************
% vedo quanti profili ci sono e leggo le variabili presenti in ciascun
% profilo

% stat_par   = ncf{'STATION_PARAMETERS'};
% s=size(stat_par);
% for i=1:s(1)
% stat_par   = ncf{'STATION_PARAMETERS'}(i,:,:)
% end
% *************************************************************************

% OPEN file
ncf = netcdf.open(file_name, 'NC_NOWRITE');

% Get all variables id
% vars = netcdf.inqVarIDs(ncf);
% for i=1:size(vars,2)
%     k = vars(i);
%     [varname, xtype, dimids, atts] = netcdf.inqVar(ncf, k);
% end
% Get variable
% [varname, xtype, dimids, atts] = netcdf.inqVar(ncid,0);
% 
% varid = netcdf.inqVarID(ncf, 'VERTICAL_SAMPLING_SCHEME');
% vsc = netcdf.getVar(ncf, varid);
% VSC   = ncf{'VERTICAL_SAMPLING_SCHEME'};inqVarId

vsc = get_var(ncf, 'VERTICAL_SAMPLING_SCHEME')';
% vsc   = get_var(ncf, 'VERTICAL_SAMPLING_SCHEME');

svsc=size(vsc);

% serve cellfun?
if cellfun(@isempty,{vsc}) || isempty(deblank(vsc))
    kv = 1;
elseif svsc(1) == 1 || svsc(2) == 1
    if strcmp('Argo CTD sampling', vsc(1,1:17)) || strcmp('Primary', vsc(1,1:7))
        kv = 1;
    else
        return
    end
else
    kv=strmatch('Primary',vsc);
    if isempty(kv)
%         keyboard
        return
    end
end
%     token = 'Primary';
%     kv = strncmp('Primary', string(vsc'), length(token));
%     ind = find(kv);
%     if isempty(ind)
%         return
%     end
% end

if svsc(1) > 1
    cmt = sprintf('MULTIPROFILE CYCLE: %s --> SAVE PRIMARY SAMPLING ONLY',  file_name(1,end-14:end));
%     disp(['MULTIPROFILE CYCLE: ' file_name(1,end-14:end) ' --> SAVE PRIMARY SAMPLING ONLY'])
    
%     app.addToLog(cmt);
    
    
end
% kv = ind(1);

if length(kv) > 1
%     keyboard
    kv=kv(1);
end


% *************************************************************************
latitude = get_var(ncf, 'LATITUDE');
latitude  = latitude(kv);

longitude = get_var(ncf, 'LONGITUDE');
longitude = longitude(kv);
 
reference_date_time = get_var(ncf, 'REFERENCE_DATE_TIME')';

dayref=datenum(sscanf(reference_date_time(1:4),'%f'),...
	    sscanf(reference_date_time(5:6),'%f'),sscanf(reference_date_time(7:8),'%f'),0,0,0);
    
juld = get_var(ncf, 'JULD');
juld = juld + dayref(kv);

%Time measurements
%======================
juld_qc = get_var(ncf, 'JULD_QC'); 
juld_qc = juld_qc(kv);

inok = [];
inok = find(double(juld_qc)~= double('1') &...
	    double(juld_qc)~= double('2') & ...
        double(juld_qc)~= double('0') & ...
        double(juld_qc)~= double('8') & ...
	    double(juld_qc)~= double('5'));
    
if ~isempty(inok)
  juld=NaN;
end

data_type   = get_var(ncf, 'DATA_TYPE')';
format_version = get_var(ncf, 'FORMAT_VERSION')';

project_name   = get_var(ncf, 'PROJECT_NAME');
pi_name  =  get_var(ncf,  'PI_NAME');


platform_number  = get_var(ncf, 'PLATFORM_NUMBER');             
cycle_number  = get_var(ncf,  'CYCLE_NUMBER');                   
direction  = get_var(ncf,  'DIRECTION');                         
data_centre  = get_var(ncf,  'DATA_CENTRE');                     
data_state_indicator  = get_var(ncf,  'DATA_STATE_INDICATOR');   
data_mode  = get_var(ncf,  'DATA_MODE');                         
dc_reference  = get_var(ncf,  'DC_REFERENCE');                   
wmo_inst_type  = get_var(ncf,  'WMO_INST_TYPE');                 

so=sort(svsc);
if so(1) > 1
    project_name = project_name(:,kv)';
    pi_name      = pi_name(:,kv)';
    platform_number = platform_number(:,kv)';
    cycle_number = cycle_number(kv);
    direction = direction(kv,:);
    data_centre = data_centre(:,kv)';
    data_state_indicator = data_state_indicator(:,kv)';
    data_mode = data_mode(:,kv)';
    dc_reference = dc_reference(:,kv)';
    wmo_inst_type = wmo_inst_type(:,kv)';
end

%Position measurements
%======================
position_qc = get_var(ncf, 'POSITION_QC');
position_qc = position_qc(kv);

inok=[];

inok =find(double(position_qc)~= double('1') &...
	   double(position_qc)~= double('2') & ...
       double(position_qc)~= double('0') & ...
       double(position_qc)~= double('8') & ...
	   double(position_qc)~= double('5'));
if ~isempty(inok)
  longitude = NaN;
  latitude = NaN;
end

%Reads measurements
%======================
pres_qc = get_var(ncf, 'PRES_QC');

if so(1) > 1
    pres_qc = pres_qc(:,kv)';
end

%PRESSURE
pressure = get_var(ncf, 'PRES');                                
if so(1) > 1   
    pressure = pressure(:,kv)';
end
pressure = qc(pressure, pres_qc);

%==================================c
pres_adjusted_qc = get_var(ncf, 'PRES_ADJUSTED_QC');  
pres_adjusted  = get_var(ncf, 'PRES_ADJUSTED');
if so(1)>1
    pres_adjusted_qc = pres_adjusted_qc(:,kv)';
    %PRESSURE_ADJUSTED
    pres_adjusted = pres_adjusted(:,kv)';
end

pres_adjusted = qc(pres_adjusted, pres_adjusted_qc);

%==================================
%TEMPERATURE
temperature  = get_var(ncf, 'TEMP');
temp_qc  = get_var(ncf, 'TEMP_QC');
if so(1)>1
    temperature = temperature(:,kv)';
    temp_qc = temp_qc(:,kv)';
end
temperature = qc(temperature, temp_qc);

%==================================
%TEMP_ADJUSTED
temp_adjusted  = get_var(ncf, 'TEMP_ADJUSTED');       
temp_adjusted_qc  = get_var(ncf, 'TEMP_ADJUSTED_QC'); 
if so(1)>1
    temp_adjusted = temp_adjusted(:,kv)';
    temp_adjusted_qc = temp_adjusted_qc(:,kv)';
end

temp_adjusted = qc(temp_adjusted, temp_adjusted_qc);

%==================================
%SALINITY
salinity  = get_var(ncf, 'PSAL');
if so(1)>1
    salinity = salinity(:,kv)';
end

if ~isempty(salinity)
    psal_qc  = get_var(ncf, 'PSAL_QC');
    if so(1)>1 
        psal_qc = psal_qc(:,kv)';
    end

    salinity = qc(salinity, psal_qc);
end

%==================================
%SAL_ADJUSTED
sal_adjusted   = get_var(ncf, 'PSAL_ADJUSTED');   
if so(1)>1
    sal_adjusted = sal_adjusted(:,kv)';
end

if ~isempty(sal_adjusted)
    psal_adjusted_qc    = get_var(ncf, 'PSAL_ADJUSTED_QC'); 
    if so(1)>1
        psal_adjusted_qc = psal_adjusted_qc(:,kv)';
    end
    
    sal_adjusted = qc(sal_adjusted, psal_adjusted_qc);
end

%==================================
% close(ncf);

clear inok;
%if T not defined => P and S not defined
inok = find(~isfinite(temperature));
if ~isempty(inok)
  pressure(inok)=NaN;
  if ~isempty(salinity)
    salinity(inok)=NaN;
  end
end

%if T not defined => P and S not defined
inok = find(~isfinite(temp_adjusted));
if ~isempty(inok)
  pres_adjusted(inok)=NaN;
  if ~isempty(sal_adjusted)
    sal_adjusted(inok)=NaN;
  end
end

iok=find(isfinite(temperature));
if isempty(iok)
  fprintf('%s\n',...
	  [' no validated (with QC=1 or QC=2 or QC=5 or QC=0 or QC=8) temperature measurements in this file']);
end

iok=find(isfinite(salinity));
if isempty(iok)
  fprintf('%s\n',...
	    [' no validated (existing or with QC=1 or QC=2 or QC=5 or QC=0 or QC=8) salinity measurements in this file ']);
end

iok=find(isfinite(temp_adjusted));
if isempty(iok)
  fprintf('%s\n',...
	  [' no validated (with QC=1 or QC=2 or QC=5 or QC=0 or QC=8) temp_adjusted measurements in this file']);
end

iok=find(isfinite(sal_adjusted));
if isempty(iok)
  fprintf('%s\n',...
	    [' no validated (existing or with QC=1 or QC=2 or QC=5 or QC=0 or QC=8) sal_adjusted measurements in this file ']);
end


% Fill output structure
argoStruct = ([]); 
argoStruct.latitude             = latitude;
argoStruct.longitude            = longitude;
argoStruct.position_qc          = position_qc;
argoStruct.juld                 = juld;
argoStruct.juld_qc              = juld_qc;
argoStruct.pressure             = pressure;
argoStruct.pres_qc              = pres_qc;
argoStruct.temperature          = temperature;
argoStruct.temp_qc              = temp_qc;
argoStruct.salinity             = salinity;
argoStruct.psal_qc              = psal_qc;
argoStruct.platform_number      = platform_number;
argoStruct.data_type            = data_type;
argoStruct.format_version       = format_version;
argoStruct.reference_date_time  = reference_date_time;
argoStruct.project_name         = project_name;
argoStruct.pi_name              = pi_name;
argoStruct.cycle_number         = cycle_number;
argoStruct.direction            = direction;
argoStruct.data_centre          = data_centre;
argoStruct.data_state_indicator = data_state_indicator;
argoStruct.data_mode            = data_mode;
argoStruct.dc_reference         = dc_reference;
argoStruct.wmo_inst_type        = wmo_inst_type;
argoStruct.pres_adjusted        = pres_adjusted;
argoStruct.pres_adjusted_qc     = pres_adjusted_qc;
argoStruct.temp_adjusted        = temp_adjusted;
argoStruct.temp_adjusted_qc     = temp_adjusted_qc;
argoStruct.sal_adjusted         = sal_adjusted;
argoStruct.psal_adjusted_qc     = psal_adjusted_qc;

% Could be a bettter way. The NetCDF file should have this info!!
[filepath,name,ext] = fileparts(file_name);

%  ind1 = strfind(file_name, 'float_data/');
%  ind2 = strfind(file_name(ind1:end), '/');
%  sl = length('float_data/');\
ind = strfind(name, '_');
argoStruct.float_name           = name(2:ind-1);
argoStruct.float_number_short   = name(ind+1:end);
argoStruct.float_number         = name;
argoStruct.file_name            = file_name;

end

% Utility functions

function data = get_var(nc, varname)
    varid = netcdf.inqVarID(nc, varname);
    data = netcdf.getVar(nc, varid);
end

function var = qc(var, var_qc)
    ind = find(double(var_qc)~= double('1') &...
	           double(var_qc)~= double('2') & ...
               double(var_qc)~= double('0') & ...
               double(var_qc)~= double('8') & ...
	           double(var_qc)~= double('5'));
    
    if ~isempty(ind)
        var(ind) = NaN;
    end
end


