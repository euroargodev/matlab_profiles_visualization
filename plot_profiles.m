function H = plot_profiles(profiles, axes)
    
    % profiles --> structure with profiles as fields
    % axes    --> axes handle to draw in

    float_num = numel(profiles);
    % =========================================================================
    hold(axes, 'on' )
    legend_names = cell(float_num,1);
    for s=1:float_num
        % get not NAN values
        jj = ~isnan(profiles{s}.salinity) & ~isnan(profiles{s}.temperature);
        H(s) = plot(axes, profiles{s}.salinity(jj), profiles{s}.temperature(jj),...
            'linewidth', 1.5, 'DisplayName', profiles{s}.float_number);
        legend_names(s,1) = {sprintf('%03d', profiles{1+float_num-s}.cycle_number)};
    end
    cmap = jet(float_num);
    cmap = flipud(cmap);
    set(H, {'color'}, num2cell(cmap, 2));

    l = legend(axes, H);
    l.Interpreter = 'none';
    l.String = legend_names;
    
    axes.XGrid = 'on';
    axes.YGrid = 'on';
    axes.XLabel.String = 'Salinity';
    axes.YLabel.String = 'Temperature (^°C)'; %'theta (^Â°C)';
    colormap(axes, 'default')
    % NOTE: modify at first show?
    header = sprintf('Float WMO %s-%03d - TS diagram, date: %s ', ...
        profiles{1}.float_name, profiles{1}.cycle_number, datestr(profiles{1}.juld(1)));
    title(axes, header);
end