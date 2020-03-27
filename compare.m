function compare(obj)
% G. Notarstefano 13 Jan 2009 - Nicola Creati 06/2019
%

basin = obj.mainApp.currentCTD;

% basin = '/data/work_giulio/giulio2/CTD_wip';

age = obj.mainApp.searchAge;
distance = obj.mainApp.searchDistance;
pressure_min = obj.mainApp.searchPressure;

% ========================================================================
% float_num = obj.currentProfile;
profile = obj.currentProfile;
float_name = profile.float_name;
float_number = profile.float_number_short;
% rootCoast = '/data/work_giulio/giulio2/COAST/';
% profile_cycle = profile.cycle_number;
pressure = profile.pressure;
salinity = profile.salinity; 
temperature = profile.temperature;

med = load(fullfile(obj.mainApp.paths.data, 'medPol.mat'));

% Convert the julian date
dateV = datestr(profile.juld(1));


% trovo le parti dei profili + profondi del valore limite (prof) per float
% e mappati% testa se il profilo è buono, occorre emettere un avviso  se
% non buono uscire
%PP=find(pressure > obj.mainApp.searchPressure);

refProfilesPath = fullfile(obj.mainApp.paths.CTD, char(basin));


% refProfilesPath = '/data/work_giulio/giulio2/CTD_wip/';

% legge i profili di riferimento
% wbox=getfname([path_2 '/*.mat']);
wbox = getfname(fullfile(refProfilesPath, '*.mat'));
[r, c] = size(wbox);
%eval(['load ' ([path_2,'map_' profile ' selected_hist'])]);

% Track all figures
figures = [];

f = figure('Name', 'Loc FP & CTDs', 'NumberTitle', 'off', 'WindowStyle', 'Docked');

medLat = med.dataset(:,1);
medLon = med.dataset(:,2);
latlim = [min(medLat(:)), max(medLat(:))];
lonlim = [min(medLon(:)), max(medLon(:))];

ax = axesm('mercator', 'MapLatLimit', latlim, 'MapLonLimit', lonlim, ...
    'Grid', 'off', 'MLineLocation', 5.0, 'MLabelLocation', 5.0,...
    'PLineLocation', 5.0, 'PLabelLocation', 5.0);
hold(ax, 'on')
geoshow(medLat, medLon);
tightmap

% Append 
figures = [figures; f];

%plot_bacini_med
spaz_time_dist=[];
sal_pre_his=[];

longitude = profile.longitude;
latitude = profile.latitude;

if longitude > 180.
    longitude = longitude -360;
end

lons = [longitude longitude+6];
lats = [latitude latitude+6];

% profile_cycle = profile.profileber_short;

timePlots = (1:5);    
timePlots(:) = NaN;      

for j=1:r

    load(fullfile(refProfilesPath, deblank(wbox(j,:))));
    % elimino profile poco profondi
    dm = size(pres);
    for kk = 1 : dm(2)
        pf = nanmax(pres(:, kk));
        if pf < pressure_min %750
            pres(:, kk) = NaN;
            sal(:, kk) = NaN;
            ptmp(:, kk) = NaN;
            temp(:, kk) = NaN;
            dates(kk) = NaN;
            lat(kk) = NaN;
            long(kk) = NaN;
            PF = 0;
            continue
        else
            PF = 1;
        end
    end
    if PF == 0
        continue
    end
    % *****************************

    mag180s = find(long > 180);
    long(mag180s) = long(mag180s) - 360;
    DATES = profile.juld(1);
    time_2 = datestr(DATES,30);
    A = strfind(time_2, 'T');
    time_2(A) = [];
    time_3 = str2double(time_2);
    DATES = changedates(time_3);
    dates_new = changedates(dates); dates_new = dates_new';
    diff_ages = dates_new-DATES; diff_age = abs(diff_ages);

    [a, b]=mercat(lons, lats);
    rng = [];
    for k = 1 : length(lat)
        range = deg2km(sqrt(((abs(long(k)) - abs(longitude))*b)^2 + ((abs(lat(k)) - abs(latitude))*b)^2));
        rng = [rng; range];
    end

    ts = find(diff_age <= age);
    ds = find(rng(ts) <= distance);

    time = dates(ts(ds))';
    latd = lat(ts(ds))';
    lond = long(ts(ds))';
    diff_t = diff_ages(ts(ds))';
    dist_s = rng(ts(ds));

    sald = sal(:,ts(ds)); 
    pred = pres(:,ts(ds));
    temd = temp(:,ts(ds));

    ss = size(sald);
    lungh=ones(ss(1),ss(2));

    for ii=1:ss(2)
        lungh(:,ii)=lungh(:,ii).*ii;
    end

    lungh = lungh(:);
    sald = sald(:);        
    pred = pred(:);
    temd = temd(:);

    lung = 1:1:length(time); lung = lung(:);

    rr = size(spaz_time_dist);
    lung = lung + rr(1);
    lungh = lungh + rr(1);

    spaz_time_dist = [spaz_time_dist; lung time lond latd diff_t dist_s];
    sal_pre_his = [sal_pre_his; lungh sald pred temd];

    % Disegna la loc dei profili in base ale differenze di tempo con
    % diversi colori
    
%     timePlots = (1:5);    
%     timePlots(:) = NaN;      
    
%     tmp = plot(ax, lond(diff_t < -9), latd(diff_t < -9),'.b', 'markersize', 14);
    tmp = geoshow(ax, latd(diff_t < -9), lond(diff_t < -9), ...
        'DisplayType','point','MarkerEdgeColor','b','Marker', '.', 'MarkerSize', 14);
    if ~isempty(tmp)
        timePlots(1) = tmp;
    end
    
%     tmp = plot(ax, lond(diff_t < -6 & diff_t >= -9),latd(diff_t < -6 & diff_t >= -9),'.c', 'markersize', 14);
    tmp = geoshow(ax, latd(diff_t < -6 & diff_t >= -9),lond(diff_t < -6 & diff_t >= -9),...
        'DisplayType','point','MarkerEdgeColor','c','Marker', '.', 'MarkerSize', 14);
    if ~isempty(tmp)
        timePlots(3) = tmp;
    end
    
 %     tmp = plot(ax, lond(diff_t < -3 & diff_t >= -6),latd(diff_t < -3 & diff_t >= -6),'.g', 'markersize', 14);
    tmp = geoshow(ax, latd(diff_t < -3 & diff_t >= -6),lond(diff_t < -3 & diff_t >= -6),...
        'DisplayType','point','MarkerEdgeColor','g','Marker', '.', 'MarkerSize', 14);
    if ~isempty(tmp)
        timePlots(2) = tmp;
    end
    
    tmp = geoshow(ax, latd(diff_t < 0 & diff_t >= -3),lond(diff_t < 0 & diff_t >= -3),...
        'DisplayType','point','MarkerEdgeColor','m','Marker', '.', 'MarkerSize', 14);
%     tmp = plot(ax, lond(diff_t < 0 & diff_t >= -3),latd(diff_t < 0 & diff_t >= -3),'.m', 'markersize', 14);
    if ~isempty(tmp)
        timePlots(4) = tmp;
    end

    tmp = geoshow(ax,  latd(diff_t >= 0),lond(diff_t >= 0),  'DisplayType','point',...
        'MarkerEdgeColor','r','Marker', '.', 'MarkerSize', 14);
%     tmp = plot(ax, lond(diff_t >= 0), latd(diff_t >= 0),'.r', 'markersize', 14 );
    if ~isempty(tmp)
        timePlots(5) = tmp;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DISEGNA LA COSTA DEL MEDITERRANEO ed il PROFILO DI RIFERIMENTO

% Disegna il float di riferimento sulla mappa del Mediterraneo
reference = geoshow(ax, latitude,longitude, 'DisplayType','point',...
    'MarkerEdgeColor','k','MarkerFaceColor','k','Marker', 'd', 'MarkerSize', 8);
% reference = plot(ax, longitude, latitude,'xk', 'markersize',8);

% Estendi timeplots per aggiunbgere la posizione del float di riferimnto
timePlots(end+1) = reference;

ddd=datestr(dateV);

tit = sprintf('Float %s - profile %s (%s) and historical CTD locations', float_name, float_number, ddd(1,1:11));
title(tit);
xlabel(ax, 'Longitude ^oE');
ylabel(ax, 'Latitude ^oN');
box(ax, 'on')

legendText = {'dt < -9 year', ...
    '-9 <= dt < -6 year',...
    '-6 <= dt < -3 year',...
    '-3 <= dt < 0 year', ...
    'time >=0 year', 'Float profile'};

ind = isnan(timePlots);
legend(ax, timePlots(~ind), legendText(~ind),'Show');
% end first plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECOND PLOT
h1=find(spaz_time_dist(:,5)>=0);
h2=find(spaz_time_dist(:,5)< 0 & spaz_time_dist(:,5) >= -3);
h3=find(spaz_time_dist(:,5)< -3 & spaz_time_dist(:,5) >= -6);
h4=find(spaz_time_dist(:,5)< -6 & spaz_time_dist(:,5) >= -9);
h5=find(spaz_time_dist(:,5)< -9);

% h=figure('parent', obj.mainFig);
% set(h, 'WindowStyle', 'Docked');

f = figure('Name', 'S FP & RP', 'NumberTitle', 'off', 'WindowStyle', 'Docked');
figures = [figures, f];
hold on

for jj=1:length(h5)
    hh=find(sal_pre_his(:,1)==h5(jj));
    %plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'.b', 'markersize',14);
    plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'-b', 'linewidth',1.2);
end 

for jj=1:length(h4)
    hh=find(sal_pre_his(:,1)==h4(jj));
    %plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'.c','markersize',14);
    plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'-c','linewidth',1.2);
end 

for jj=1:length(h3)
    hh=find(sal_pre_his(:,1)==h3(jj));
    %plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'.g', 'markersize',14);
    plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'-g', 'linewidth',1.2);
end 

for jj=1:length(h2)
    hh=find(sal_pre_his(:,1)==h2(jj));
    %plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'.m', 'markersize',14);
    plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'-m', 'linewidth',1.2);
end 

for jj=1:length(h1)
    hh=find(sal_pre_his(:,1)==h1(jj));
    %plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'.r', 'markersize',14);
    plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'-r', 'linewidth',1.2);
end

%ref = plot(salinity,-pressure,'.k', 'markersize',14);
ref = plot(salinity,-pressure,'-k', 'linewidth',1.8);

titolo = sprintf('Sal profiles: float %s - profile: %s (%s) vs historical ctd', float_name, float_number, ddd(1,1:11));
% title(['Sal profiles: float # ' profile.float_name ' - profile ' num2str(profile) ' (' ddd(1,1:11) ') vs historical ctd']);
title(titolo);
xlabel('Salinity');
ylabel('Pressure (dbar)');

timePlots(end) = ref;

ms=nanmin(salinity(:)); Ms=nanmax(salinity(:));
Mp=nanmax(pressure(:));

set(gca,'xlim',[roundn(ms,-2)-0.25 roundn(Ms,-2)+0.25],'ylim',[-(roundn(Mp,+2))-100 0]); % ,'xtick',[37.1:0.2:39.5] ,'ytick',[-2500:200:0]
box on
grid on

legendText = {'dt < -9 year', ...
    '-9 <= dt < -6 year',...
    '-6 <= dt < -3 year',...
    '-3 <= dt < 0 year', ...
    'time >=0 year', 'Float profile'};

ind = isnan(timePlots);

legend(gca, timePlots(~ind), legendText(~ind),'Show');
% end second plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SECOND PLOT TS
% h1=find(spaz_time_dist(:,5)>=0);
% h2=find(spaz_time_dist(:,5)< 0 & spaz_time_dist(:,5) >= -3);
% h3=find(spaz_time_dist(:,5)< -3 & spaz_time_dist(:,5) >= -6);
% h4=find(spaz_time_dist(:,5)< -6 & spaz_time_dist(:,5) >= -9);
% h5=find(spaz_time_dist(:,5)< -9);

% h=figure('parent', obj.mainFig);
% set(h, 'WindowStyle', 'Docked');

f = figure('Name', 'TS FP & RP', 'NumberTitle', 'off', 'WindowStyle', 'Docked');
figures = [figures, f];
hold on

for jj=1:length(h5)
    hh=find(sal_pre_his(:,1)==h5(jj));
    %plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'.b', 'markersize',14);
    plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'-b', 'linewidth',1.2);
end 

for jj=1:length(h4)
    hh=find(sal_pre_his(:,1)==h4(jj));
    %plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'.c','markersize',14);
    plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'-c', 'linewidth',1.2);
end 

for jj=1:length(h3)
    hh=find(sal_pre_his(:,1)==h3(jj));
    %plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'.g', 'markersize',14);
    plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'-g', 'linewidth',1.2);
end 

for jj=1:length(h2)
    hh=find(sal_pre_his(:,1)==h2(jj));
    %plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'.m', 'markersize',14);
    plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'-m', 'linewidth',1.2);
end 

for jj=1:length(h1)
    hh=find(sal_pre_his(:,1)==h1(jj));
    %plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'.r', 'markersize',14);
    plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'-r', 'linewidth',1.2);
end

%ref = plot(salinity,temperature,'.k', 'markersize',14);
ref = plot(salinity,temperature,'-k', 'linewidth',1.8);

titolo = sprintf('TS profiles: float %s - profile: %s (%s) vs historical ctd', float_name, float_number, ddd(1,1:11));
% title(['Sal profiles: float # ' profile.float_name ' - profile ' num2str(profile) ' (' ddd(1,1:11) ') vs historical ctd']);
title(titolo);
xlabel('Salinity (PSU)');
ylabel('Temperature (^oC)');

timePlots(end) = ref;

ms=nanmin(salinity(:)); Ms=nanmax(salinity(:));
mp=nanmin(temperature(:)); Mp=nanmax(temperature(:));

set(gca,'xlim',[roundn(ms,-2)-0.5 roundn(Ms,-2)+0.25],'ylim',[roundn(mp,-2)-1 roundn(Mp,-2)+1]); % ,'xtick',[37.1:0.2:39.5] ,'ytick',[-2500:200:0]
box on
grid on

legendText = {'dt < -9 year', ...
    '-9 <= dt < -6 year',...
    '-6 <= dt < -3 year',...
    '-3 <= dt < 0 year', ...
    'time >=0 year', 'Float profile'};

ind = isnan(timePlots);

legend(gca, timePlots(~ind), legendText(~ind),'Show');
% end second plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THIRD PLOT

[mi,im]=sort(abs(spaz_time_dist(:,5)));

for giu=1:length(im)
    hh=find(sal_pre_his(:,1)==im(giu));
    nd=isnan(sal_pre_his(hh,2));
    if all(nd)
        continue
    else
        im=im(giu);
        break
    end
end

timediff=spaz_time_dist(im,5);  % frazione di anno (devo moltiplicare * 365 il valore ottenuto)

% FIGURA PROFILO FLOAT CON PROFILO PIU' VICINO IN TEMPO
f = figure('Name', 'Salinity Profiles dt', 'NumberTitle', 'off', 'WindowStyle', 'Docked');
figures = [figures, f];
hold on
box on
grid on

hh=find(sal_pre_his(:,1)==im);
%plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'.r', 'markersize', 14);
plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'-r', 'linewidth', 1.2);

AAA=spaz_time_dist(im,2);
AA=num2str(AAA);
A=datevec(AA(1:8),'yyyymmdd');
DA=datenum([num2str(A(3)),'-',num2str(A(2)),'-',num2str(A(1))],'dd-mm-yyyy');
DAA=datestr(DA);

%plot(salinity, -pressure,'.k', 'markersize', 14);
plot(salinity, -pressure,'-k', 'linewidth', 1.8);

titolo = sprintf('Float %s - profile: %s (%s) vs nearest (in time) ctd (%s)', float_name, float_number, ddd(1,1:11), DAA);
title(titolo);

% title(['Float ' profile ' - profile ' num2str(profile) ' (' ddd(1,1:11) ') vs nearest (in time) ctd (' DAA ')']);
xlabel('Salinity (PSU)');
ylabel('Pressure (dbar)');

ms=nanmin(salinity(:)); Ms=nanmax(salinity(:));
Mp=nanmax(pressure(:));

set(gca,'xlim',[roundn(ms,-2)-0.25 roundn(Ms,-2)+0.25],'ylim',[-(roundn(Mp,+2))-100 0]); % ,'xtick',[37.1:0.2:39.5] ,'ytick',[-2500:200:0]

%set(gca,'xlim',[38.4 39.6],'xtick',[38.4:0.2:39.6],'ylim',[-2500 0],'ytick',[-2500:250:0]);

l = sprintf('time diff (hist-float) = %s days', num2str(round(timediff*365)));
dim = [.5 .5 .3 .3];
annotation('textbox',dim,'String',l,'FitBoxToText','on');
% end thid plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% THIRD PLOT TS

% [mi,im]=sort(abs(spaz_time_dist(:,5)));
% 
% for giu=1:length(im)
%     hh=find(sal_pre_his(:,1)==im(giu));
%     nd=isnan(sal_pre_his(hh,2));
%     if all(nd)
%         continue
%     else
%         im=im(giu);
%         break
%     end
% end
% 
% timediff=spaz_time_dist(im,5);  % frazione di anno (devo moltiplicare * 365 il valore ottenuto)

% FIGURA PROFILO FLOAT CON PROFILO PIU' VICINO IN TEMPO
f = figure('Name', 'TS Profiles dt', 'NumberTitle', 'off', 'WindowStyle', 'Docked');
figures = [figures, f];
hold on
box on
grid on

hh=find(sal_pre_his(:,1)==im);
%plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'.r', 'markersize', 14);
plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'-r', 'linewidth', 1.2);

AAA=spaz_time_dist(im,2);
AA=num2str(AAA);
A=datevec(AA(1:8),'yyyymmdd');
DA=datenum([num2str(A(3)),'-',num2str(A(2)),'-',num2str(A(1))],'dd-mm-yyyy');
DAA=datestr(DA);

plot(salinity, temperature,'-k', 'linewidth', 1.8);

titolo = sprintf('Float %s - profile: %s (%s) vs nearest (in time) ctd (%s)', float_name, float_number, ddd(1,1:11), DAA);
title(titolo);

% title(['Float ' profile ' - profile ' num2str(profile) ' (' ddd(1,1:11) ') vs nearest (in time) ctd (' DAA ')']);
xlabel('Salinity (PSU)');
ylabel('Temperature (^oC)');

ms=nanmin(salinity(:)); Ms=nanmax(salinity(:));
mp=nanmin(temperature(:)); Mp=nanmax(temperature(:));

set(gca,'xlim',[roundn(ms,-2)-0.5 roundn(Ms,-2)+0.25],'ylim',[roundn(mp,-2)-1 roundn(Mp,-2)+1]); % ,'xtick',[37.1:0.2:39.5] ,'ytick',[-2500:200:0

%set(gca,'xlim',[roundn(ms,-2)-0.25 roundn(Ms,-2)+0.25],'ylim',[-(roundn(Mp,+2))-100 0]); % ,'xtick',[37.1:0.2:39.5] ,'ytick',[-2500:200:0]

%set(gca,'xlim',[38.4 39.6],'xtick',[38.4:0.2:39.6],'ylim',[-2500 0],'ytick',[-2500:250:0]);

l = sprintf('time diff (hist-float) = %s days', num2str(round(timediff*365)));
dim = [.5 .5 .3 .3];
annotation('textbox',dim,'String',l,'FitBoxToText','on');
% end thid plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FOURTH PLOT
f = figure('Name', 'MAP dt ', 'NumberTitle', 'off', 'WindowStyle', 'Docked');
figures = [figures, f];
ax = axesm('mercator', 'MapLatLimit', latlim, 'MapLonLimit', lonlim, ...
    'Grid', 'on', 'MLineLocation', 5.0, 'MLabelLocation', 5.0,...
    'PLineLocation', 5.0, 'PLabelLocation', 5.0);
hold(ax, 'on')
geoshow(medLat, medLon)
tightmap

geoshow(spaz_time_dist(im,4),spaz_time_dist(im,3), 'DisplayType','point','MarkerEdgeColor','r','Marker', 'd', 'MarkerSize', 8);
geoshow(latitude, longitude,'DisplayType','point','MarkerEdgeColor','k','Marker', 'd', 'MarkerSize', 8);

titolo = sprintf('Float %s - profile: %s (%s) and nearest (in time) ctd (%s)', float_name, float_number, ddd(1,1:11), DAA);
title(titolo);
xlabel('Longitude ^oE');
ylabel('Latitude ^oN');

spacediff=spaz_time_dist(im,6);
l = sprintf('distance = %s km', num2str(round(spacediff)));
dim = [.2 0.01 .3 .3];
annotation('textbox', dim, 'String', l, 'FitBoxToText','on');

[mi,im]=sort(abs(spaz_time_dist(:,6)));

for giu=1:length(im)
    hh=find(sal_pre_his(:,1)==im(giu));
    nd=isnan(sal_pre_his(hh,2));
    if all(nd)
        continue
    else
        im=im(giu);
        break
    end
end
% end fourth plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIFTH PLOT
% FIGURA PROFILO FLOAT CON PROFILO PIU' VICINO IN SPAZIO
timediff=spaz_time_dist(im,5);  % frazione di anno (devo moltiplicare * 365 il valore ottenuto)
%[mi,im]=min(abs(spaz_time_dist(:,6)));

f = figure('Name', 'Salinity Profiles ds ', 'NumberTitle', 'off', 'WindowStyle', 'Docked');
figures = [figures, f];
hold on

hh=find(sal_pre_his(:,1)==im);
%plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'.r', 'markersize',14);
plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'-r', 'linewidth', 1.2);

AAA=spaz_time_dist(im,2);
AA=num2str(AAA);
A=datevec(AA(1:8),'yyyymmdd');
DA=datenum([num2str(A(3)),'-',num2str(A(2)),'-',num2str(A(1))],'dd-mm-yyyy');
DAA=datestr(DA);

%plot(salinity,-pressure,'.k', 'markersize',14);
plot(salinity,-pressure,'-k', 'linewidth', 1.8);

titolo = sprintf('Float %s - profile: %s (%s) vs nearest (in space) ctd (%s)', float_name, float_number, ddd(1,1:11), DAA);
title(titolo);

xlabel('Salinity (PSU)');
ylabel('Pressure (dbar)');

ms=nanmin(salinity(:)); Ms=nanmax(salinity(:));
Mp=nanmax(pressure(:));

set(gca,'xlim',[roundn(ms,-2)-0.25 roundn(Ms,-2)+0.25],'ylim',[-(roundn(Mp,+2))-100 0]); % ,'xtick',[37.1:0.2:39.5] ,'ytick',[-2500:200:0]

%set(gca,'xlim',[38.4 39.6],'xtick',[38.4:0.2:39.6],'ylim',[-2500 0],'ytick',[-2500:250:0]);
box on
grid on

l = sprintf('time diff (hist.-float) = %s days', num2str(round(timediff*365)));
dim = [.2 .5 .3 .3];
annotation('textbox', dim, 'String', l, 'FitBoxToText','on');
% end fifth plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FIFTH PLOT TS
% FIGURA PROFILO FLOAT CON PROFILO PIU' VICINO IN SPAZIO
%[mi,im]=min(abs(spaz_time_dist(:,6)));

f = figure('Name', 'TS Profiles ds ', 'NumberTitle', 'off', 'WindowStyle', 'Docked');
figures = [figures, f];
hold on

hh=find(sal_pre_his(:,1)==im);
%plot(sal_pre_his(hh,2),-sal_pre_his(hh,3),'.r', 'markersize',14);
plot(sal_pre_his(hh,2),sal_pre_his(hh,4),'-r', 'linewidth', 1.2);

AAA=spaz_time_dist(im,2);
AA=num2str(AAA);
A=datevec(AA(1:8),'yyyymmdd');
DA=datenum([num2str(A(3)),'-',num2str(A(2)),'-',num2str(A(1))],'dd-mm-yyyy');
DAA=datestr(DA);

%plot(salinity,-pressure,'.k', 'markersize',14);
plot(salinity,temperature,'-k', 'linewidth', 1.8);

titolo = sprintf('Float %s - profile: %s (%s) vs nearest (in space) ctd (%s)', float_name, float_number, ddd(1,1:11), DAA);
title(titolo);

xlabel('Salinity (PSU)');
ylabel('Temperature (^oC)');

ms=nanmin(salinity(:)); Ms=nanmax(salinity(:));
mp=nanmin(temperature(:)); Mp=nanmax(temperature(:));

set(gca,'xlim',[roundn(ms,-2)-0.5 roundn(Ms,-2)+0.25],'ylim',[roundn(mp,-2)-1 roundn(Mp,-2)+1]); % ,'xtick',[37.1:0.2:39.5] ,'ytick',[-2500:200:0

%set(gca,'xlim',[roundn(ms,-2)-0.25 roundn(Ms,-2)+0.25],'ylim',[-(roundn(Mp,+2))-100 0]); % ,'xtick',[37.1:0.2:39.5] ,'ytick',[-2500:200:0]

%set(gca,'xlim',[38.4 39.6],'xtick',[38.4:0.2:39.6],'ylim',[-2500 0],'ytick',[-2500:250:0]);
box on
grid on

l = sprintf('time diff (hist.-float) = %s days', num2str(round(timediff*365)));
dim = [.2 .5 .3 .3];
annotation('textbox', dim, 'String', l, 'FitBoxToText','on');
% end fifth plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIXTH PLOT
f = figure('Name', 'MAP ds ', 'NumberTitle', 'off', 'WindowStyle', 'Docked');
figures = [figures, f];
hold on
box on
ax = axesm('mercator', 'MapLatLimit', latlim, 'MapLonLimit', lonlim, ...
    'Grid', 'on', 'MLineLocation', 5.0, 'MLabelLocation', 5.0,...
    'PLineLocation', 5.0, 'PLabelLocation', 5.0);
hold(ax, 'on')
geoshow(medLat, medLon)
tightmap

geoshow(spaz_time_dist(im,4),spaz_time_dist(im,3), 'DisplayType','point','MarkerEdgeColor','r','Marker', 'd', 'MarkerSize', 8);
geoshow(latitude, longitude ,'DisplayType','point','MarkerEdgeColor','k','Marker', 'd', 'MarkerSize', 8);

titolo = sprintf('Float %s - profile: %s  (%s) and nearest (in space) ctd (%s)', float_name, float_number, ddd(1,1:11), DAA);
title(titolo);
xlabel('Longitude ^oE');
ylabel('Latitude ^oN');

spacediff=spaz_time_dist(im,6);

l = sprintf('distance  = %s km', num2str(round(spacediff)));
dim = [.2 0.1 .3 .3];
annotation('textbox', dim, 'String', l, 'FitBoxToText','on');

obj.mainApp.figures = figures;

% time difference output
% timediff*365

