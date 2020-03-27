% I path vanno decisi ed impostati in qualche modo !!!
classdef ArgoUI < matlab.apps.AppBase
   properties
      btn_Compare
      btn_loadArgoData
      btn_pickArgoData
      btn_plotArgoData 
      btn_Reset
      comboRef
      currentCTD
      figures
      log
      logText
      mainFig
      onPicking
      paths
      plot
      pn_controls
      profiles
      profilesInfo
      profilesList
      refProfile  
      searchAge
      searchDistance
      searchPressure
      spn_Age
      spn_Distance      
      spn_Pressure  
      winPosition
   end
   
    methods
           % Set app paths
        function setPaths(obj, paths)
            obj.paths = paths;
        end

        % Class constructor
        function obj = ArgoUI(appPaths)
           obj.mainFig = uifigure('Name','Argo', 'Color', 'none', ...
            'Visible', 'on','Position', [0 0 430 620],'Resize', 'off');%, ...
            %'CloseRequestFcn',@(fig, event) obj.onClose(fig)); 
            
            % workaround to avoid wrong starting positioning and sizing of
            % the UI
            %set(0,'defaultfigureposition',[100 100 100 100]);
            
            obj.winPosition = obj.mainFig.Position;
            obj.onPicking = 0;
            obj.paths = appPaths;
            obj.currentCTD = string('adriatico');
            obj.searchDistance = 250e3;
            obj.searchPressure = 250e3;
            obj.searchAge = 10;
            obj.setUpControl() 
            obj.figures = [];
        end
        
        % Setup every control
        function setUpControl(obj)
            pn_position = obj.mainFig.Position;
            [x0, y0, x1, y1] = feval(@(x) x{:}, num2cell([obj.mainFig.Position]));
            % Load float data
            obj.btn_loadArgoData = uibutton('parent', obj.mainFig, ...
                'Text', 'Load data','Position', [x0+10, obj.winPosition(4)-30, 200, 22],...
                'ButtonPushedFcn', @(btn,event) obj.onLoadArgoData(btn));
            % Plot floats
            obj.btn_plotArgoData = uibutton('parent', obj.mainFig, 'Text', 'Plot data', ...
                'Position', [x0+10, obj.winPosition(4)-60, 200, 22],...
                'Enable', 'off','ButtonPushedFcn', @(btn,event) obj.onPlotArgoData(btn));  
            % Pick float profiles
            obj.btn_pickArgoData = uibutton(obj.mainFig, 'state', 'Text', 'Pick curve', ...
                'Position', [x0+10, obj.winPosition(4)-90, 200, 22],...
                'Enable', 'off', 'ValueChangedFcn', @(btn,event) obj.onPickingArgoPlot(btn));              
            % Set the reference profile
            obj.comboRef = uidropdown('parent', obj.mainFig, 'Editable', 'off', ...
                'Enable', 'off','Position', [x0+10, obj.winPosition(4)-120, 200, 22],...
                'ValueChangedFcn', @(comboRef, event )obj.selectCompareDataset(comboRef));
            % Select age parameter
            obj.spn_Age = uispinner(obj.mainFig, 'Value', 10, 'Step', 10, ...
                'ValueDisplayFormat','%u year', 'Limits', [1, 100], ...
                'Position', [x0+10, obj.winPosition(4)-150, 200, 22],...
                'Enable', 'off',...
                'ValueChangedFcn', @(spn_Age, event) obj.onSelectAge(spn_Age));     
            % Select distance (in km) parameter
            obj.spn_Distance = uispinner(obj.mainFig, 'Value', 250, ...
                'Step', 10, 'ValueDisplayFormat','%u km','Limits', [1, 1000], ...
                'Position', [x0+10, obj.winPosition(4)-180, 200, 22],...
                'Enable', 'off',...
                'ValueChangedFcn', @(spn_Distance, event) obj.onSelectDistance(spn_Distance));         
            % Select pressure (in dbar) parameter
            obj.spn_Pressure = uispinner(obj.mainFig, 'Value', 750, ...
                'Step', 100, 'ValueDisplayFormat','%u dbar','Limits', [0, 1000], ...
                'Position', [x0+10, obj.winPosition(4)-210, 200, 22],...
                'Enable', 'off',...
                'ValueChangedFcn', @(spn_Pressure, event) obj.onSelectPressure(spn_Pressure));             
            % Compare reference to current selected profile
            obj.btn_Compare = uibutton(obj.mainFig, 'Text', 'Compare', ...
                'Position', [x0+10, obj.winPosition(4)-240, 200, 22],...
                'Enable', 'off','ButtonPushedFcn', @(btn,event) obj.onCompare(btn));  
            % Show current selected float info
            obj.profilesInfo = uitextarea('parent', obj.mainFig,'Value', 'Info',...
                 'Editable', 'off', 'Position', [x0+10, obj.winPosition(4)-410, 200, 150]);  
            % Show loaded floats
            obj.profilesList = uilistbox('parent', obj.mainFig, 'Items',{'Profiles'},...
                 'Position', [x0+220, obj.winPosition(4)-410, 200, 402], ...
                 'ValueChangedFcn', @(lbx,event) obj.onSelectionChanged(lbx));             
            % Log area                         
            obj.log = uitextarea('parent', obj.mainFig,'Editable', 'off', ...
                 'Position', [x0+10, obj.winPosition(4)-590, 410, 170]);  
            obj.log.Value='Log';

            obj.btn_Reset = uibutton('parent', obj.mainFig, 'Text', 'Reset', ...
                'Position', [x0+10, obj.winPosition(4)-615, 410, 22],...
                'Enable', 'on', 'BackgroundColor', [1.0, 0.0, 0.0],...
                'FontColor', [1 1 1],...
                'ButtonPushedFcn', @(btn,event) obj.onReset(btn));  
            
            % Fill compare dropBox
            ctd = dir(obj.paths.CTD);

            items = cell(1,length(ctd)-2);
            for i=3:length(ctd)   
                items(1,i-2) = {ctd(i).name};   
            end
            obj.comboRef.Items = items;
                                     
        end      
        
        function onReset(obj, evt)
            if ~isempty(obj.figures)
                delete(obj.figures)
            end
            obj.onPicking = 0;            
            obj.currentCTD = string();
            obj.figures = []; 
            obj.plot.enablePicking(0)
            set(obj.btn_plotArgoData, 'Enable', 'off');
            delete(obj.plot.mainFig);
            obj.profiles = struct();
            
            obj.profilesList.set('Items', {});  
            % Enable controls
            set(obj.btn_pickArgoData, 'Enable', 'off');
            set(obj.btn_Compare, 'Enable', 'off');
            set(obj.spn_Distance, 'Enable', 'off');
            set(obj.spn_Pressure, 'Enable', 'off');
            set(obj.spn_Age, 'Enable', 'off');
            set(obj.comboRef, 'Enable', 'off'); 
            set(obj.profilesInfo,'Value', '');
            
            
        end
        
        function onSelectDistance(obj, evt)
            obj.searchDistance = evt.Value;
        end

        function onSelectPressure(obj, evt)
            obj.searchPressure = evt.Value;
        end        
        
        function onSelectAge(obj, evt)
            obj.searchAge = evt.Value;
        end
        
        
        function selectCompareDataset(obj, event)
            obj.currentCTD = event.Value;
        end
        
        function setDefaultProfiles(obj, names)
            obj.refProfile.Items = names;
        end
        
        function onSelectionChanged(obj, lbx)            
            indexCell = strfind(lbx.Items, lbx.Value);
            index = find(not(cellfun('isempty', indexCell)));
            H = obj.plot.axes.Children;
            obj.plot.selectProfiles(H(index)) ;
        end
               
        function onClose(obj, fig)
            % Close all figures
            selection = questdlg('Close Argo?',...
                'Confirmation',...
                'Yes','No','Yes');
            switch selection,
                case 'Yes',
                    delete(get(obj.mainFig, 'Children'));
                    delete(obj.mainFig);
%                     delete(get(groot, 'Children'));
%                     set(0,'ShowHiddenHandles','on')
%                     delete(get(0,'Children'))

                case 'No'
                    return
            end
        end
        
        function onPickingArgoPlot(obj, btn, evt)
            obj.plot.enablePicking(btn.Value)
        end
    
        
        function onCompare(obj, btn, evt)
            obj.plot.compareProfile()
        end        
    
        
        function onLoadArgoData(obj, btn)
            disp('Leggo i dati')
            obj.readData()
            set(obj.btn_plotArgoData, 'Enable', 'on');
%             set(obj.btn_pickArgoData, 'Enable', 'on');
        end

        function onPlotArgoData(obj, btn)
            obj.plot = ArgoPlot(obj);

            a= cell(1,length(obj.profiles));
            for i=1:length(obj.profiles)
              a(i)= {obj.profiles{i}.float_number};
            end
            obj.profilesList.set('Items', a);  
            % Enable controls
            set(obj.btn_pickArgoData, 'Enable', 'on');
            set(obj.btn_Compare, 'Enable', 'on');
            set(obj.spn_Distance, 'Enable', 'on');
            set(obj.spn_Pressure, 'Enable', 'on');
            set(obj.spn_Age, 'Enable', 'on');
            set(obj.comboRef, 'Enable', 'on');
        end        
        
        function readData(obj)
            % Permette di decidere i file da filtrare
            tt = {'*.nc', 'NetCDF File'; '*.mat', 'MAT-Files'};

            % List of profiles to read and precess
            pfs = uipickfiles('FilterSpec', obj.paths.float_data, 'Type', tt, 'Prompt', 'Select Profile');


            % Profile data cell container
            obj.profiles = cell(length(pfs), 1);

            % Read process every selected profiles
            for i=1:length(obj.profiles)
                obj.profiles{i} = process_profiles(pfs{i}, obj);
            end
            
        end
        
        function addToLog(obj, str)
            currString = obj.logText.Value;
            currString{end+1}=str;
%             set(obj.logText, 'Value', currString);
            obj.logText.Value = currString;
        end
            

    end
end