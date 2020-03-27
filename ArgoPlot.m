% I path vanno decisi ed impostati in qualche modo !!!
% Devo usare la funzione Figure standard datro che i nuovi oggetti ui non
% catturano alcune eventi!!!
classdef ArgoPlot < handle
   properties
      mainApp
      mainFig
      profiles
      axes
%       winPosition
%       toolbar
      currentProfile
%       mainFig2
   end
   
    methods           
        % Class constructor
        function obj = ArgoPlot(app)
            obj.mainApp = app;
            obj.mainFig = figure('Name','TS', ...
                       'Visible', 'on', ... 
                       'Color', 'white', ...
                       'WindowStyle', 'docked', ...
                       'Numbertitle', 'off',...
                       'CloseRequestFcn',@(fig, event) obj.onClose(fig));
%             
%            obj.compareFig = figure('Name','Compare plot', ...                        
%                        'Visible', 'off', ...
%                        'WindowStyle', 'docked', ...
%                        'Numbertitle', 'off',...
%                        'WindowButtonDownFcn', @obj.mapPick);
           
%             obj.winPosition = obj.mainFig.Position;
            obj.profiles = app.profiles;

            obj.axes = axes('Parent', obj.mainFig);
            plot_profiles(obj.profiles, obj.axes); 
            
            % Get the figure toolbar handle
            obj.currentProfile = NaN;

        end 
        
        function mapPick(obj, fig, evt)
                mousePos=get(fig.CurrentAxes,'CurrentPoint');
                disp(['You clicked X:',num2str(mousePos(1)),',  Y:',num2str(mousePos(2))]);
        end
        
        function onClose(obj, fig)
            fig.Visible = 'off';
        end     
       
        function compareProfile(obj)
            compare(obj)
            
        end        
        
        function enablePicking(obj, value)
            lines = obj.axes.Children;
            for i=1:length(lines)
                if value
                    set(lines(i), 'ButtonDownFcn', @obj.MouseLineSelected)   
                else
                    set(lines(i), 'ButtonDownFcn', {})   
                end
            end
        end
                 
        function MouseLineSelected(obj, line, event)
            obj.selectProfiles(line)
        end

        function selectProfiles(obj, line)
             set(line, 'LineWidth', 2.5);
             H = obj.axes.Children;
             set(H(H ~= line), 'LineWidth', 0.5);
             ind = find(H==line);
             
             obj.currentProfile = obj.profiles{ind};
             
             listBox = obj.mainApp.profilesList;             
             
             listBox.Value = obj.currentProfile.float_number;
             
             % Change the title to the right profile
             header = sprintf('Float WMO %s-%03d - TS diagram, date: %s, ', ...
             obj.currentProfile.float_name, obj.currentProfile.cycle_number, ...
             datestr(datestr(obj.currentProfile.juld(1))));
             title(obj.axes, header);             
            
             %cellArrayText = cell(1,5);
             cellArrayText = cell(1,6);
             cellArrayText{1} = sprintf('Project Name: \t%s', obj.currentProfile.project_name);
             cellArrayText{2} = sprintf('Float Name: \t%s', obj.currentProfile.float_name);
             cellArrayText{3} = sprintf('Cycle: \t\t%d', obj.currentProfile.cycle_number);
             cellArrayText{4} = sprintf('Longitude: \t%f', obj.currentProfile.longitude);
             cellArrayText{5} = sprintf('Latitude: \t%f', obj.currentProfile.latitude);
             cellArrayText{6} = sprintf('Time: \t%s', datestr(obj.currentProfile.juld(1)));
             obj.mainApp.profilesInfo.set('Value', cellArrayText);


        end   
    end
end
