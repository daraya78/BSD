classdef uiBSD < handle
    
    properties
        gui
        study
        default_eeg
        is_saved
        ws_filename
    end
    
    methods
        function obj=uiBSD()
            % loaded data
            obj.study.data = {};
            obj.study.metadata = {};
            % generated data
            obj.study.generated_data = {};
            obj.study.generated_metadata = {};
            % viterbi sequence
            obj.study.viterbi_results = {};
            obj.study.viterbi_metadata = {};
            % training results
            obj.study.training_results = {};
            obj.study.training_metadata = {};
            % decoded results
            obj.study.decode_results = {};
            obj.study.decode_metadata = {};
            % model
            obj.study.model = {};
            % current study filename
            obj.ws_filename = '';
            obj.is_saved = 1;
            
            obj.default_eeg = {'ele_bio64.xyz', 'ele_bio128.xyz'};
            
            % ui creation
            obj.gui.fig = figure(...
                'units', 'pixels',...
                'position', [100 100 600 400],...
                'Toolbar', 'none',...
                'Menubar', 'none',...
                'name', 'BSD',...
                'numbertitle', 'off',...
                'resize', 'off');

            obj.update_title();
            % Menues
            obj.create_file_menu();
            obj.create_model_menu();
            obj.create_inference_menu();
            obj.create_view_menu();
            
            obj.create_main_screen();
                        
            % show study status
            obj.show_study();
            
            obj.load_opt();
            
            obj.gui.fig.CloseRequestFcn = @obj.leave;
        end
        
        function load_opt(obj)
            program_option;
            obj.study.opt = opt;
        end
        
        % create menues
        function create_file_menu(obj)
            % File menu
            obj.gui.file_menu = uimenu(obj.gui.fig,...
                'Label', 'Data');
            
            % Import Data submenu
            obj.gui.import_menu = uimenu(obj.gui.file_menu,...
                'Label', 'Import data');
            % Open Ascii Action
            obj.gui.open_ascii_menu = uimenu(obj.gui.import_menu,...
                'Label', 'From ASCII file');
            % Open Mat File Action
            obj.gui.open_matlab_menu = uimenu(obj.gui.import_menu,...
                'Label', 'From Matlab file');
            % Open Global Workspace Variable Action
            obj.gui.open_matlab_ws_menu = uimenu(obj.gui.import_menu,...
                'Label', 'From Matlab Workspace');

            % Export Data submenu
            obj.gui.export_menu = uimenu(obj.gui.file_menu,...
                'Label', 'Export generated data',...
                'Enable', 'off');
            % Save Ascii Action
            obj.gui.save_ascii_menu = uimenu(obj.gui.export_menu,...
                'Label', 'To ASCII file');
            % Save Mat File Action
            obj.gui.save_matlab_menu = uimenu(obj.gui.export_menu,...
                'Label', 'To Matlab file');
            % Save to Matlab workspace Action
            obj.gui.save_matlab_ws_menu = uimenu(obj.gui.export_menu,...
                'Label', 'To Matlab Workspace');

            % Create new study complete
            obj.gui.new_study_menu = uimenu(obj.gui.file_menu,...
                'Label', 'Create New Study',...
                'Separator', 'on');
            % Load study complete
            obj.gui.load_study_complete_menu = uimenu(obj.gui.file_menu,...
                'Label', 'Load Study',...
                'Separator', 'off');
            % Load study
            obj.gui.load_study_partial_menu = uimenu(obj.gui.file_menu,...
                'Label', 'Load Items from Study',...
                'Separator', 'off');
            % Save study
            obj.gui.save_study_menu = uimenu(obj.gui.file_menu,...
                'Label', 'Save Current Study',...
                'Enable', 'off');
            % Save as study
            obj.gui.save_as_study_menu = uimenu(obj.gui.file_menu,...
                'Label', 'Save Current Study As',...
                'Enable', 'off');
            % clear study
            obj.gui.clear_study_menu = uimenu(obj.gui.file_menu,...
                'Label', 'Clear Study',...
                'Enable', 'off');
            % quit app
            obj.gui.quit_menu = uimenu(obj.gui.file_menu,...
                'Label', 'Quit',...
                'Separator', 'on');

            % set callback
            % import
            obj.gui.open_ascii_menu.Callback = @obj.open_ascii_action;
            obj.gui.open_matlab_menu.Callback = @obj.open_matlab_action;
            obj.gui.open_matlab_ws_menu.Callback = @obj.open_matlab_ws_action;
            % export
            obj.gui.save_ascii_menu.Callback = @obj.save_ascii_action;
            obj.gui.save_matlab_menu.Callback = @obj.save_matlab_action;
            obj.gui.save_matlab_ws_menu.Callback = @obj.save_matlab_ws_action;
            % other menues
            obj.gui.new_study_menu.Callback = @obj.new_study_action;
            obj.gui.load_study_complete_menu.Callback = {@obj.load_study_action, 0};
            obj.gui.load_study_partial_menu.Callback = {@obj.load_study_action, 1};
            obj.gui.save_study_menu.Callback = {@obj.save_study_action, 0};
            obj.gui.save_as_study_menu.Callback = {@obj.save_study_action, 1};
            obj.gui.clear_study_menu.Callback = @obj.clear_study_action;
            obj.gui.quit_menu.Callback = @obj.leave;
        end
           
        function create_model_menu(obj)
            % Model menu
            obj.gui.model_menu = uimenu(obj.gui.fig,...
                'Label', 'Model');
            % Create Model Action
            obj.gui.create_model_menu = uimenu(obj.gui.model_menu,...
                'Label', 'Create Model');
            % Prior Distribution Action
            obj.gui.prior_dist_menu = uimenu(obj.gui.model_menu,...
                'Label', 'Prior Distributions',...
                'Enable', 'off');
            % Posteriori Distribution Action
            obj.gui.posterior_dist_menu = uimenu(obj.gui.model_menu,...
                'Label', 'Posterior Distributions',...
                'Enable', 'off');
            
            % set callbacks
            obj.gui.create_model_menu.Callback = @obj.create_model_action;
            obj.gui.prior_dist_menu.Callback = @obj.prior_dist_action;
            obj.gui.posterior_dist_menu.Callback = @obj.posteriori_dist_action;
        end

        function create_inference_menu(obj)
            % Inference menu
            obj.gui.inference_menu = uimenu(obj.gui.fig,...
                'Label', 'Inference');
            % Training Model Action
            obj.gui.train_model_menu = uimenu(obj.gui.inference_menu,...
                'Label', 'Train Model',...
                'Enable', 'off');
            % Estimate Model Parameters Action
            obj.gui.estimate_parameters_menu = uimenu(obj.gui.inference_menu,...
                'Label', 'Estimate Parameters',...
                'Enable', 'off');
            % Data Decode Action
            obj.gui.data_decode_menu = uimenu(obj.gui.inference_menu,...
                'Label', 'State Probabilities',...
                'Enable', 'off');
            % Generate from Model Action
            obj.gui.generate_data_menu = uimenu(obj.gui.inference_menu,...
                'Label', 'Generate Data',...
                'Enable', 'off');
            % Sequence estimation by Viterbi
            obj.gui.sequence_estimation_menu = uimenu(obj.gui.inference_menu,...
                'Label', 'Viterbi Decoding',...
                'Enable', 'off');
            
            % set callbacks
            obj.gui.train_model_menu.Callback = @obj.train_model_action;
            obj.gui.estimate_parameters_menu.Callback = @obj.estimate_parameters_action;
            obj.gui.data_decode_menu.Callback = @obj.data_decode_action;
            obj.gui.generate_data_menu.Callback = @obj.generate_data_action;
            obj.gui.sequence_estimation_menu.Callback = @obj.sequence_estimation_action;
        end

        function create_view_menu(obj)
            % View menu
            obj.gui.view_menu = uimenu(obj.gui.fig,...
                'Label', 'View');
            % Secuence Action
            obj.gui.plot_sequence_menu = uimenu(obj.gui.view_menu,...
                'Label', 'States Sequence',...
                'Enable', 'off');
            % Ocupancia(poner en ingles) Action
            obj.gui.plot_occupancy_menu = uimenu(obj.gui.view_menu,...
                'Label', 'Occupancy',...
                'Enable', 'off');
            % duration histogram Action
            obj.gui.plot_histogram_menu = uimenu(obj.gui.view_menu,...
                'Label', 'Duration Histogram',...
                'Enable', 'off');
            % duration histogram Action
            obj.gui.plot_topo_menu = uimenu(obj.gui.view_menu,...
                'Label', 'Topography',...
                'Enable', 'off');
            % duration histogram Action
            obj.gui.plot_cov_matrix_menu = uimenu(obj.gui.view_menu,...
                'Label', 'Emission Covariance Matrix',...
                'Enable', 'off');
            % duration histogram Action
            obj.gui.plot_free_energy_menu = uimenu(obj.gui.view_menu,...
                'Label', 'Free Energy',...
                'Enable', 'off');
            % gamma plot
            obj.gui.plot_gamma_menu = uimenu(obj.gui.view_menu,...
                'Label', 'State Probability',...
                'Enable', 'off');
            % gamma plot
            obj.gui.plot_posterior_dur_menu = uimenu(obj.gui.view_menu,...
                'Label', 'Duration Posterior',...
                'Enable', 'off');
            
            % generated data view menu
            obj.gui.gen_view_menu = uimenu(obj.gui.view_menu,...
                'Label', 'Generated data');
            % generated data view menu
            obj.gui.plot_gen_sequence_menu = uimenu(obj.gui.gen_view_menu,...
                'Label', 'States Sequence',...
                'Enable', 'off');
            % generated data view menu
            obj.gui.plot_gen_occupancy_menu = uimenu(obj.gui.gen_view_menu,...
                'Label', 'Occupancy',...
                'Enable', 'off');
            % generated data view menu
            obj.gui.plot_gen_durhist_menu = uimenu(obj.gui.gen_view_menu,...
                'Label', 'Duration Histogram',...
                'Enable', 'off');
            % generated data view menu
            obj.gui.plot_gen_topo_menu = uimenu(obj.gui.gen_view_menu,...
                'Label', 'Topography',...
                'Enable', 'off');
            % generated data view menu
            obj.gui.plot_gen_covmat_menu = uimenu(obj.gui.gen_view_menu,...
                'Label', 'Emmision Covariance Matrix',...
                'Enable', 'off');
            
            % set callbacks
            obj.gui.plot_sequence_menu.Callback = @obj.plot_sequence_action;
            obj.gui.plot_occupancy_menu.Callback = @obj.plot_occupancy_action;
            obj.gui.plot_histogram_menu.Callback = @obj.plot_histogram_action;
            obj.gui.plot_topo_menu.Callback = @obj.plot_topo_action;
            obj.gui.plot_cov_matrix_menu.Callback = @obj.plot_cov_matrix_action;
            obj.gui.plot_free_energy_menu.Callback = @obj.plot_free_energy;
            obj.gui.plot_gamma_menu.Callback = @obj.plot_gamma;
            obj.gui.plot_posterior_dur_menu.Callback = @obj.plot_posterior_dur;
            % generated data callbacks
            obj.gui.plot_gen_sequence_menu.Callback = @obj.plot_gen_sequence_action;
            obj.gui.plot_gen_occupancy_menu.Callback = @obj.plot_gen_occupancy_action;
            obj.gui.plot_gen_durhist_menu.Callback = @obj.plot_gen_histogram_action;
            obj.gui.plot_gen_topo_menu.Callback = @obj.plot_gen_topo_action;
            obj.gui.plot_gen_covmat_menu.Callback = @obj.plot_gen_cov_matrix_action;
        end

        % create main screen
        function create_main_screen(obj)
            % tab container
            obj.gui.tabs_container = uitabgroup(obj.gui.fig);

            % tab for status
            obj.gui.tab_status = uitab(obj.gui.tabs_container, 'Title', 'Status');

            % create status bar
            obj.gui.status_bar = uicontrol(obj.gui.fig,...
                'units', 'pixels',...
                'Position', [10 5 558 25],...
                'FontSize', 12,...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', 'Ready');

            y = 10;
            % temp
            y = y + 165;
            obj.gui.model_panel = uipanel(obj.gui.tab_status,...
                'TitlePosition', 'centertop',...
                'units', 'pixels',...
                'Position', [10 y 578 105],...
                'FontSize', 10,...
                'Title', 'Model');
            obj.gui.prior_label = uicontrol(obj.gui.model_panel,...
                'units', 'pixels',...
                'Position', [10 60 135 15],...
                'FontSize', 9,...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', 'Prior Distributions');
            obj.gui.prior_value = uicontrol(obj.gui.model_panel,...
                'units', 'pixels',...
                'Position', [145 60 90 15],...
                'FontSize', 9,...
                'FontWeight', 'bold',...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', '???');
            obj.gui.posterior_label = uicontrol(obj.gui.model_panel,...
                'units', 'pixels',...
                'Position', [10 35 135 15],...
                'FontSize', 9,...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', 'Posterior Distributions');
            obj.gui.posterior_value = uicontrol(obj.gui.model_panel,...
                'units', 'pixels',...
                'Position', [145 35 90 15],...
                'FontSize', 9,...
                'FontWeight', 'bold',...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', '???');
            obj.gui.nstates_label = uicontrol(obj.gui.model_panel,...
                'units', 'pixels',...
                'Position', [10 10 60 15],...
                'FontSize', 9,...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', '# States');
            obj.gui.nstates_value = uicontrol(obj.gui.model_panel,...
                'units', 'pixels',...
                'Position', [70 10 90 15],...
                'FontSize', 9,...
                'FontWeight', 'bold',...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', '???');
            
            y = y + obj.gui.model_panel.Position(4);
            obj.gui.data_panel = uipanel(obj.gui.tab_status,...
                'TitlePosition', 'centertop',...
                'units', 'pixels',...
                'Position', [10 y 578 80],...
                'FontSize', 10,...
                'Title', 'Data');
            % filename label and edit
            obj.gui.data_filename_label = uicontrol(obj.gui.data_panel,...
                'units', 'pixels',...
                'Position', [10 35 60 15],...
                'FontSize', 9,...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', 'Data File:');
            obj.gui.data_filename_value = uicontrol(obj.gui.data_panel,...
                'units', 'pixels',...
                'Position', [70 35 458 15],...
                'FontSize', 9,...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', '');
            % dimension label and edit
            obj.gui.data_dimension_label = uicontrol(obj.gui.data_panel,...
                'units', 'pixels',...
                'Position', [10 10 90 15],...
                'FontSize', 9,...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', '# of Channels:');
            obj.gui.data_dimension_value = uicontrol(obj.gui.data_panel,...
                'units', 'pixels',...
                'Position', [100 10 40 15],...
                'FontSize', 9,...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', '');
            % length label and edit
            obj.gui.data_length_label = uicontrol(obj.gui.data_panel,...
                'units', 'pixels',...
                'Position', [180 10 75 15],...
                'FontSize', 9,...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', 'Data Length:');
            obj.gui.data_length_value = uicontrol(obj.gui.data_panel,...
                'units', 'pixels',...
                'Position', [265 10 70 15],...
                'FontSize', 9,...
                'HorizontalAlignment', 'left',...
                'Style', 'text',...
                'String', '');
        end

        % callbacks
        function open_ascii_action(obj, ~, ~)
            % ask for file
            [selected_file, path] = uigetfile({'*.*',  'All Files (*.*)'});
            % if a file is selected
            if selected_file ~= 0
                % open and store it
                filename = [path selected_file];
                % read data
                try
                    data = load(filename);
                catch err
                    obj.my_error(err);
                    return;
                end
                % store in the object
                obj.study = setfield(obj.study, 'data', data);
                % ask for a name
                name = inputdlg({'Name the data:'}, 'Name?', [1 50], {'my data'});
                % if cancel or give empty string assing default one
                if isempty(name)
                    name = {'my data'};
                end
                % store metadata
                metadata = {}
                metadata.filename = [path selected_file];
                metadata.name = name{1};
                metadata.dimension = size(data, 2);
                metadata.length = size(obj.study.data, 1);
                metadata.multisubject = 0;
                obj.study = setfield(obj.study, 'metadata', metadata);
                % refresh info 
                obj.show_study();
                obj.set_status_bar(['Loaded "' selected_file '"']);
            end
            obj.set_saved_status(0);
        end

        function open_matlab_action(obj, ~, ~)
            % ask for file
            [selected_file, path] = uigetfile({'*.mat',  'Matlab File (*.mat)'});
            % if a file is selected
            if selected_file ~= 0
                % open the data
                filename = [path selected_file];
                try
                    new_data = load(filename);
                catch err
                    obj.my_error(err);
                    return;
                end
                % get the keys
                fields_data = fields(new_data);
                % display a selection list
                [indx,tf] = listdlg('PromptString', 'Select Data Variable:',...
                    'ListSize', [250 150],...
                    'SelectionMode', 'single',...
                    'ListString', fields_data);
                if ~tf
                    % cancel the loading by close or press cancel
                    return;
                end
                % get the first field
                value = getfield(new_data, fields_data{indx});
                
                % store metadata
                metadata = {};
                % ask if data is a struct (multisubject)
                metadata.multisubject = 0;
                if isstruct(value)
                    metadata.multisubject = 1;
                end
                
                % store the data
                obj.study.data = value;
                % ask for a name
                name = inputdlg({'Name the data:'}, 'Name?', [1 50], {'my data'});
                % if cancel or give empty string assing default one
                if isempty(name)
                    name = {'my data'};
                end
                
                metadata.filename = strcat(path, selected_file);
                metadata.name = name{1};
                if metadata.multisubject
                    metadata.dimension = size(obj.study.data.cond(1).subj(1).block(1).data, 2);
                    metadata.length = {};
                else
                    metadata.dimension = size(obj.study.data, 2);
                    metadata.length = size(obj.study.data, 1);
                end
                obj.study.metadata = metadata;
                % refresh info 
                obj.show_study();
                obj.set_status_bar(['Loaded "' fields_data{indx} '" from file "' selected_file '" Loaded']);
            end
            obj.set_saved_status(0);
        end

        function open_matlab_ws_action(obj, ~, ~)
            % retrieve var from global study
            vars = evalin('base', 'who()');
            % ask the user to select one
            [indx,tf] = listdlg('PromptString', 'Select Variable Data:',...
                'ListSize', [250 150],...
                'SelectionMode', 'single',...
                'ListString', vars);
            if ~tf
                % cancel the loading by close or press cancel
                return;
            end
            % get the value
            value = evalin('base', vars{indx});
            
            % store metadata
            metadata = {};
            % ask if data is a struct (multisubject)
            metadata.multisubject = 0;
            if isstruct(value)
                metadata.multisubject = 1;
            end

            % assign data
            obj.study.data = value;
            metadata.filename = '(from matlab workspace)';
            metadata.name = vars{indx};
            if metadata.multisubject
                metadata.dimension = size(obj.study.data.cond(1).subj(1).block(1).data, 2);
                metadata.length = {};
            else
                metadata.dimension = size(obj.study.data, 2);
                metadata.length = size(obj.study.data, 1);
            end
            obj.study.metadata = metadata;
            % refresh info 
            obj.show_study();
            obj.set_status_bar(['Loaded "' vars{indx} '" from Workspace']);
            obj.set_saved_status(0);
        end

        function save_ascii_action(obj, ~, ~)
            % ask filename
            [selected_file, path] = uiputfile({'*.*',  'All Files (*.*)'});
            filename = [path selected_file];
            try
                % define format
                format = [repmat('%.6f ', 1, size(obj.study.generated_data, 2)) '\n'];
                % open file
                fileID = fopen(filename, 'w');
                for i=1:size(obj.study.generated_data,1)
                    % write every line
                    fprintf(fileID, format, obj.study.generated_data(i,:));
                end
                % close file
                fclose(fileID);
                obj.set_status_bar(['Created "' selected_file '"']);
            catch err
                obj.my_error(err);
                return;
            end
        end
        
        function save_matlab_action(obj, ~, ~)
            % ask filename
            [selected_file, path] = uiputfile({'*.mat',  'Matlab File (*.mat)'});
            filename = [path selected_file];
            % auxiliar variable
            data = obj.study.generated_data;
            %store filename in matlab file
            try
                save(filename, 'data');
                obj.set_status_bar(['Created "' selected_file '"']);
            catch err
                obj.my_error(err);
                return;
            end
        end
        
        function save_matlab_ws_action(obj, ~, ~)
            name = inputdlg({'Name variable:'}, 'Name?', [1 50], {'generated_data'});
            if isempty(name)
                name = {'generated_data'};
            end
            assignin('base', name{1}, obj.study.generated_data);
            obj.set_status_bar(['Set "' name{1} ' in the Workspace']);
        end
        
        function new_study_action(obj, ~, ~)
            % first ask for ask the previous study
            if ~obj.is_save_query('Save?', 'Save previous Study?')
                return;
            end
            % loaded data
            obj.study.data = {};
            obj.study.metadata = {};
            % generated data
            obj.study.generated_data = {};
            obj.study.generated_metadata = {};
            % viterbi sequence
            obj.study.viterbi_results = {};
            obj.study.viterbi_metadata = {};
            % training results
            obj.study.training_results = {};
            obj.study.training_metadata = {};
            % decoded results
            obj.study.decode_results = {};
            obj.study.decode_metadata = {};
            % model
            obj.study.model = {};
            % current study filename
            obj.ws_filename = '';
            % load again opt
            obj.load_opt();
            % update all
            obj.set_status_bar(['New Study created']);
            obj.set_saved_status(0);
        end
        
        function result=save_study_action(obj, ~, ~, save_as)
            result=0;
            data = obj.get_select_vars('Save');
            if (~isempty(data))
                % if at least one of the study fields is selected
                something_to_save = any(cell2mat(arrayfun(@(t){getfield(data, t{1});}, fields(data))));
                if ~something_to_save
                    obj.my_warning('Must select something');
                    result=0;
                    return;
                end
                
                % ask for file if save_as
                if save_as || isempty(obj.ws_filename)
                    [selected_file, path] = uiputfile({'*.mat',  'Matlab File (*.mat)'});
                    % if a file is selected
                    if selected_file == 0
                        result=0;
                        return;
                    end
                    
                    obj.ws_filename = [path selected_file];
                end

                % opt values
                if isfield(data, 'opt') && data.opt
                    ws.opt = obj.study.opt;
                else
                    obj.study.opt = {};
                end
                % viterbi sequence
                if isfield(data, 'viterbi') && data.viterbi
                    ws.viterbi_results = obj.study.viterbi_results;
                    ws.viterbi_metadata = obj.study.viterbi_metadata;
                else
                    obj.study.viterbi_results = {};
                    obj.study.viterbi_metadata = {};
                end
                % generated data
                if isfield(data, 'generated') && data.generated
                    ws.generated_data = obj.study.generated_data;
                    ws.generated_metadata = obj.study.generated_metadata;
                else
                    obj.study.generated_data = {};
                    obj.study.generated_metadata = {};
                end
                % decoded data
                if isfield(data, 'decode') && data.decode
                    ws.decode_results = obj.study.decode_results;
                    ws.decode_metadata = obj.study.decode_metadata;
                else
                    obj.study.decode_results = {};
                    obj.study.decode_metadata = {};
                end
                % results
                if isfield(data, 'results') && data.results
                    ws.training_results = obj.study.training_results;
                    ws.training_metadata = obj.study.training_metadata;
                else
                    obj.study.training_results = {};
                    obj.study.training_metadata = {};
                end
                % model
                if isfield(data, 'model') && data.model
                    ws.model = obj.study.model;
                else
                    obj.study.model = {};
                end
                % training data
                if isfield(data, 'data') && data.data
                    ws.data = obj.study.data;
                    ws.metadata = obj.study.metadata;
                else
                    obj.study.data = {};
                    obj.study.metadata = {};
                end

                try
                    save(obj.ws_filename, 'ws');
                catch err
                    obj.my_error(err);
                    result=0;
                    return;
                end
                pieces = strsplit(obj.ws_filename, '\');
                obj.set_status_bar(['Study "' pieces{size(pieces,2)} '" Saved']);
                obj.set_saved_status(1);
                result=1;
            end
        end

        function load_study_action(obj, ~, ~, partial)
            % first ask for ask the previous study
            if ~obj.is_save_query('Save?', 'Save previous Study?')
                return;
            end
            
            % ask for file to load
            [selected_file, path] = uigetfile({'*.mat',  'Matlab File (*.mat)'});
            % if a file is selected
            if selected_file ~= 0
                % open the data
                filename = [path selected_file];
                try
                    data = load(filename);
                catch err
                    obj.my_error(err);
                    return;
                end
                % drop old stuff if not partial
                if ~partial
                    % loaded data
                    obj.study.data = {};
                    obj.study.metadata = {};
                    % generated data
                    obj.study.generated_data = {};
                    obj.study.generated_metadata = {};
                    % viterbi sequence
                    obj.study.viterbi_results = {};
                    obj.study.viterbi_metadata = {};
                    % training results
                    obj.study.training_results = {};
                    obj.study.training_metadata = {};
                    % decoded results
                    obj.study.decode_results = {};
                    obj.study.decode_metadata = {};
                    % model
                    obj.study.model = {};
                end
                    
                % new one is the loaded one
                new_path = filename;

                % put variables in study
                try
                    vars = fields(data.ws);
                catch err
                    obj.my_error(err);
                    return;
                end
                v = {};
                if partial
                    v = obj.get_select_vars2('Load', vars);
                else
                    for i=1:size(vars,1)
                        v = setfield(v, vars{i}, 1);
                    end
                end
                
                % if nothing selected to load just leave
                something_to_load = any(cell2mat(arrayfun(@(t){getfield(v, t{1});}, fields(v))));
                if ~something_to_load
                    if partial
                        obj.my_warning('Must select something');
                    else
                        obj.my_warning('Empty Study');
                    end
                    return;
                end
                % options
                if isfield(v, 'opt') && v.opt
                    var_value = getfield(data.ws, 'opt');
                    obj.study = setfield(obj.study, 'opt', var_value);
                end
                % viterbi
                if (isfield(v, 'viterbi') && v.viterbi) || (isfield(v, 'viterbi_results') && v.viterbi_results)
                    var_value = getfield(data.ws, 'viterbi_results');
                    obj.study = setfield(obj.study, 'viterbi_results', var_value);
                    var_value = getfield(data.ws, 'viterbi_metadata');
                    obj.study = setfield(obj.study, 'viterbi_metadata', var_value);
                end
                % generated_data
                if (isfield(v, 'generated') && v.generated) || (isfield(v, 'generated_data') && v.generated_data)
                    var_value = getfield(data.ws, 'generated_data');
                    obj.study = setfield(obj.study, 'generated_data', var_value);
                    var_value = getfield(data.ws, 'generated_metadata');
                    obj.study = setfield(obj.study, 'generated_metadata', var_value);
                end
                % decode_results
                if (isfield(v, 'decode') && v.decode) || (isfield(v, 'decode_results') && v.decode_results)
                    var_value = getfield(data.ws, 'decode_results');
                    obj.study = setfield(obj.study, 'decode_results', var_value);
                    var_value = getfield(data.ws, 'decode_metadata');
                    obj.study = setfield(obj.study, 'decode_metadata', var_value);
                end
                % training_results
                if (isfield(v, 'results') && v.results) || (isfield(v, 'training_results') && v.training_results)
                    var_value = getfield(data.ws, 'training_results');
                    obj.study = setfield(obj.study, 'training_results', var_value);
                    var_value = getfield(data.ws, 'training_metadata');
                    obj.study = setfield(obj.study, 'training_metadata', var_value);
                end
                % model
                if isfield(v, 'model') && v.model
                    var_value = getfield(data.ws, 'model');
                    obj.study = setfield(obj.study, 'model', var_value);
                end
                % data
                if isfield(v, 'data') && v.data
                    var_value = getfield(data.ws, 'data');
                    obj.study = setfield(obj.study, 'data', var_value);
                    var_value = getfield(data.ws, 'metadata');
                    obj.study = setfield(obj.study, 'metadata', var_value);
                end
                
                obj.ws_filename = filename;
                obj.show_study();
                obj.set_status_bar(['Study "' selected_file '" Loaded']);
                obj.set_saved_status(1);
            end
        end

        function clear_study_action(obj, ~, ~)
            data = obj.get_select_vars('Clear');
            if (~isempty(data))
                % opt values
                if isfield(data, 'opt') && data.opt
                    obj.study.opt = {};
                end
                % viterbi sequence
                if isfield(data, 'viterbi') && data.viterbi
                    obj.study.viterbi_results = {};
                    obj.study.viterbi_metadata = {};
                end
                % generated data
                if isfield(data, 'generated') && data.generated
                    obj.study.generated_data = {};
                    obj.study.generated_metadata = {};
                end
                % decoded data
                if isfield(data, 'decode') && data.decode
                    obj.study.decode_results = {};
                    obj.study.decode_metadata = {};
                end
                % results
                if isfield(data, 'results') && data.results
                    obj.study.training_results = {};
                    obj.study.training_metadata = {};
                end
                % model
                if isfield(data, 'model') && data.model
                    obj.study.model = {};
                end
                % training data
                if isfield(data, 'data') && data.data
                    obj.study.data = {};
                    obj.study.metadata = {};
                end
                % drop the filename
                obj.ws_filename = '';
                
                obj.show_study();
                obj.set_status_bar('Study Dropped');
                obj.set_saved_status(0);
            end
        end
        
        function create_model_action(obj, ~, ~)
            if ~obj.check_opt()
                return;
            end
            
            data = obj.get_model();
            drawnow();

            if (~isempty(data))
                obj.study.model = data.model;
                obj.show_study();
                obj.set_status_bar('Model Created');
                % reset generated data
                obj.study.generated_data = {};
                obj.study.generated_metadata = {};
                % reset viterbi sequence
                obj.study.viterbi_results = {};
                obj.study.viterbi_metadata = {};
                % reset training results
                obj.study.training_results = {};
                obj.study.training_metadata = {};
                % reset decoded results
                obj.study.decode_results = {};
                obj.study.decode_metadata = {};
                % set as unsaved
                obj.set_saved_status(0);
            end
        end
        
        function prior_dist_action(obj, ~, ~)
            data = obj.get_distribution(1);
            drawnow();
            
            if (~isempty(data))
                % put the new parameters in the model
                if data.multisubject
                    model_name = class(obj.study.model.matrixmodel(1,1,1))
                    obj.study.model.matrixmodel(1,1,1).emis_model.prior = data.model.cond(1).subj(1).block(1).emis;
                    if strcmp(model_name, 'hsmm')
                        obj.study.model.matrixmodel(1,1,1).dur_model.prior = data.model.cond(1).subj(1).block(1).dur;
                    end
                    obj.study.model.matrixmodel(1,1,1).trans_model.prior = data.model.cond(1).subj(1).block(1).trans;
                    obj.study.model.matrixmodel(1,1,1).in_model.prior = data.model.cond(1).subj(1).block(1).init;
                else
                    model_name = class(obj.study.model);
                    obj.study.model.emis_model.prior = data.model.emis;
                    if strcmp(model_name, 'hsmm')
                        obj.study.model.dur_model.prior = data.model.dur;
                    end
                    obj.study.model.trans_model.prior = data.model.trans;
                    obj.study.model.in_model.prior = data.model.init;
                end
                obj.set_status_bar('Prior Modified');
                obj.set_saved_status(0);
            end
        end
           
        function posteriori_dist_action(obj, ~, ~)
            obj.get_distribution(0);
            drawnow();
        end

        function train_model_action(obj, ~, ~)
            if ~obj.check_opt()
                return;
            end

            data = obj.get_estimation();
            drawnow();
            
            if (~isempty(data))
                if isempty(obj.study.model)
                    obj.my_warning('No model created'); 
                    return
                end
                if isempty(obj.study.data)
                    obj.my_warning('No data loaded'); 
                    return
                end
                if strcmp(class(obj.study.model), 'multimodel')
                    obj.my_warning('Multimodel cannot train, please create a new model');
                    return
                end
                
                obj.set_status_bar('Training Model...');
                % set learning algorithm
                if strcmp(data.learning_algorithm, 'Variational Bayes')
                    obj.study.opt = setfield(obj.study.opt, 'train', 'VB');
                elseif strcmp(data.learning_algorithm, 'Expectation Maximization ')
                    obj.study.opt = setfield(obj.study.opt, 'train', 'EM');
                end
                % set init option
                if strcmp(data.initialization, 'K-means')
                    obj.study.opt = setfield(obj.study.opt, 'initoption', 'kmeans');
                elseif strcmp(data.initialization, 'Random')
                    obj.study.opt = setfield(obj.study.opt, 'initoption', 'random');
                end
                
                % set numeric parameters
                obj.study.opt = setfield(obj.study.opt, 'matixer', data.maxiter);
                obj.study.opt = setfield(obj.study.opt, 'nrep', data.n_repetitions);
                obj.study.opt = setfield(obj.study.opt, 'tol', data.tolerance);
                obj.study.opt = setfield(obj.study.opt, 'maxitersim', data.max_iter_sim);
                obj.study.opt = setfield(obj.study.opt, 'maxcyc', data.max_cycles);
                obj.study.opt = setfield(obj.study.opt, 'dmax', data.dmax);
                obj.study.opt = setfield(obj.study.opt, 'minstates', data.n_states{1});
                obj.study.opt = setfield(obj.study.opt, 'maxstates', data.n_states{2});
                
                if data.multisubject
                    obj.study.opt = setfield(obj.study.opt, 'shareemis', data.shareemis);
                    obj.study.opt = setfield(obj.study.opt, 'sharein', data.sharein);
                    obj.study.opt = setfield(obj.study.opt, 'sharetrans', data.sharetrans);
                    obj.study.opt = setfield(obj.study.opt, 'sharedur', data.sharedur);
                end
                    
                % set init option
                if data.verbose
                    obj.study.opt = setfield(obj.study.opt, 'verbose', 'yes');
                else
                    obj.study.opt = setfield(obj.study.opt, 'verbose', 'no');
                end
                try
                    if data.multisubject
                        [results, devnull1, devnull2, new_model] = obj.study.model.train(obj.study.data, 'opt', obj.study.opt);
                        obj.study.model = new_model;
                        obj.study.training_results = results;
                    else
                        results = obj.study.model.train(obj.study.data, 'opt', obj.study.opt);
                        obj.study.training_results = results;
                    end
                    obj.study.training_metadata.multisubject = data.multisubject;
                    obj.show_study();
                    obj.set_status_bar('Model Trained');
                    % reset generated data
                    obj.study.generated_data = {};
                    obj.study.generated_metadata = {};
                    % reset viterbi sequence
                    obj.study.viterbi_results = {};
                    obj.study.viterbi_metadata = {};
                    % reset decoded results
                    obj.study.decode_results = {};
                    obj.study.decode_metadata = {};
                    % set as unsaved
                    obj.set_saved_status(0);
                catch err
                    obj.my_error(err);
                    obj.set_status_bar('Model Training Failed');
                end
            end
        end

        function estimate_parameters_action(obj, ~, ~)
            % ask if exit both options
            t_r = ~isempty(obj.study.training_results);
            v_r = ~isempty(obj.study.viterbi_results);
            
            if t_r && v_r
                choice = questdlg('Do you want to use model data or viterbi data?',...
                    'Estimate Parameters',...
                    'model', 'viterbi', 'model');
                if isempty(choice)
                    t_r = 0;
                    v_r = 0;
                end
                if strcmp(choice, 'model')
                    t_r = 1;
                    v_r = 0;
                end
                if strcmp(choice, 'viterbi')
                    t_r = 0;
                    v_r = 1;
                end
                
                % if nothing exist or nothing selected just return 
                if ~t_r && ~v_r
                    return;
                end
                
                if t_r
                    obj.study.model.emis_model.update(obj.study.opt.train, obj.study.data, obj.study.training_results.decodevar.gamma);
                    obj.study.model.trans_model.update(obj.study.opt.train, obj.study.training_results.decodevar.xi);
                    obj.study.model.in_model.update(obj.study.opt.train, obj.study.training_results.decodevar.gamma(1,:))
                    obj.study.model.dur_model.update(obj.study.opt.train, (1:opt.dmax)', obj.study.training_results.decodevar.durcount);
                else
                    obj.study.model.emis_model.update(obj.study.opt.train, obj.study.data, obj.study.viterbi_results.decodevar.gamma);
                    obj.study.model.trans_model.update(obj.study.opt.train, obj.study.viterbi_results.decodevar.xi);
                    obj.study.model.in_model.update(obj.study.opt.train, obj.study.viterbi_results.decodevar.gamma(1,:))
                    obj.study.model.dur_model.update(obj.study.opt.train, (1:opt.dmax)', obj.study.viterbi_results.decodevar.durcount);
                end
                obj.set_status_bar('Parameters Estimated');
                obj.set_saved_status(0);
            end
        end

        function data_decode_action(obj, ~, ~)
            if ~obj.check_opt()
                return;
            end
            
            choice = questdlg('Do you want to calculate state probabilities?',...
                'State Probabilities',...
                'Yes', 'No', 'Yes');
            if strcmp(choice, 'Yes')
                obj.set_status_bar('State probabilities in progress...');
                try
                    decode_results = obj.study.model.decode(obj.study.data, obj.study.opt);
                    obj.study.decode_results = decode_results;
                    obj.study.decode_metadata.multisubject = obj.study.metadata.multisubject;
                    obj.set_status_bar('State probabilities done');
                    obj.set_saved_status(0);
                catch err
                    obj.my_error(err);
                    obj.set_status_bar('State probabilities Failed');
                end
            end
        end
        
        function generate_data_action(obj, ~, ~)
            if ~obj.check_opt()
                return;
            end
            
            if strcmp(class(obj.study.model), 'multimodel')
                obj.my_warning('What should the program do when you "generate" data with a multisubject model?');
                return;
            end
            
            if isempty(obj.study.training_results) && obj.study.model.nstates == 0
                obj.my_warning('Model must have a defined number of states to generate data');
                return;
            end
            
            data = obj.get_generation();
            drawnow();
            
            if (~isempty(data))
                obj.set_status_bar('Generating data...');
                model_name = class(obj.study.model);
                % set data
                obj.study.model.trans_model.parsampl = data.parameters.trans;
                obj.study.model.in_model.parsampl = data.parameters.in;
                obj.study.model.emis_model.parsampl = data.parameters.emis;
                if strcmp(model_name, 'hsmm')
                    obj.study.model.dur_model.parsampl = data.parameters.dur;
                end
                % generate data
                try
                    [new_data seq] = obj.study.model.gen(data.n);
                    obj.study.generated_data = new_data;
                    obj.study.generated_metadata.seq = seq;
                    obj.set_status_bar('Data generated');
                    obj.set_saved_status(0);
                catch err
                    obj.my_error(err);
                    obj.set_status_bar('Data generation Failed');
                end
            end
        end
        
        function sequence_estimation_action(obj, ~, ~)
            if ~obj.check_opt()
                return;
            end

            choice = questdlg('Do you want to get most probable sequence?',...
                'Viterbi Decoding',...
                'Yes', 'No', 'Yes');
            if strcmp(choice, 'Yes')
                obj.set_status_bar('Calculating Most Probable Sequence...');
                try
                    viterbi_results = obj.study.model.viterbi(obj.study.data, obj.study.opt);
                    obj.study.viterbi_results = viterbi_results;
                    obj.study.viterbi_metadata.multisubject = obj.study.metadata.multisubject;
                    obj.set_status_bar('Most Probable Sequence Calculated');
                    obj.set_saved_status(0);
                catch err
                    obj.my_error(err);
                    obj.set_status_bar('Most Probable Sequence Failed');
                end
            end
        end

        function plot_sequence_action(obj, ~, ~)
            selected = [1, 1, 1];

            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'States Sequence',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            function my_draw(c, s, b)
                if strcmp(class(obj.study.model), 'multimodel')
                    my_values = obj.study.viterbi_results.cond(c).subj(s).block(b).stateseq;
                else
                    my_values = obj.study.viterbi_results;
                end
                b = bar(my_values);
                xlabel('Time');
                ylabel('State');
                title('State Label');
            end
            
            try
                my_draw(selected(1), selected(2), selected(3));
                if obj.study.viterbi_metadata.multisubject
                    obj.multisubject_selector(f, @my_draw);
                    set(f.CurrentAxes, 'Units', 'pixels');
                    set(f.CurrentAxes, 'Position', [f.CurrentAxes.Position(1) f.CurrentAxes.Position(2) f.CurrentAxes.Position(3) f.CurrentAxes.Position(4)-50]);
                end
            catch err
                obj.my_error(err);
            end
        end
        
        function plot_occupancy_action(obj, ~, ~)
            selected = [1, 1, 1];
            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'Occupancy',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            function my_draw(c, s, b)
                if strcmp(class(obj.study.model), 'multimodel')
                    my_model = obj.study.model.matrixmodel(c, s, b);
                    seq = obj.study.viterbi_results.cond(c).subj(s).block(b).stateseq;
                else
                    my_model = obj.study.model;
                    seq = obj.study.viterbi_results;
                end
                [a b]=hist(seq,1:my_model.nstates);
                bar(b,100.0*a/size(seq,2));
                xlabel('State Label');
                ylabel('% Ocupancy');
                title('State Ocupancy');

            end
            
            try
                my_draw(selected(1), selected(2), selected(3));
                if obj.study.viterbi_metadata.multisubject
                    obj.multisubject_selector(f, @my_draw);
                    set(f.CurrentAxes, 'Units', 'pixels');
                    set(f.CurrentAxes, 'Position', [f.CurrentAxes.Position(1) f.CurrentAxes.Position(2) f.CurrentAxes.Position(3) f.CurrentAxes.Position(4)-50]);
                end
            catch err
                obj.my_error(err);
            end
        end
        
        function plot_histogram_action(obj, ~, ~)
            selected = [1, 1, 1];
            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'Duration Histograms',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            function my_draw(c, s, b)
                if strcmp(class(obj.study.model), 'multimodel')
                    my_model = obj.study.model.matrixmodel(c, s, b);
                    seq = obj.study.viterbi_results.cond(c).subj(s).block(b).stateseq;
                else
                    my_model = obj.study.model;
                    seq = obj.study.viterbi_results;
                end
                
                [a b]=util.calculadur(seq);
                for j=1:my_model.nstates
                    subplot(my_model.nstates,1,j);
                    [s1 s2]=hist(a(b==j));
                    hist(a(b==j),obj.study.opt.dmax);
                    axis([0 obj.study.opt.dmax 0 max(s1)+1])
                    ylabel('Counts');
                    xlabel('State Duration');
                    %Agregado David

                    if j==1
                        title('State Duration Histograms');
                    end
                end
            end
            
            try
                my_draw(selected(1), selected(2), selected(3));
                if obj.study.viterbi_metadata.multisubject
                    obj.multisubject_selector(f, @my_draw);
                end
            catch err
                obj.my_error(err);
            end
        end

        function plot_topo_action(obj, ~, ~)
            is_tp = which('topoplot');
            if isempty(is_tp)
                obj.my_warning('It requires Eeglab');
                delete(f);
                return;
            end
            
            ruta = obj.select_eeg_coords();
            if isempty(ruta)
                obj.my_warning('Must choose a coordinates file for your eeg hardware');
                return;
            end
            
            selected = [1, 1, 1];
            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'Topography',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            function my_draw(c, s, b)
                if strcmp(class(obj.study.model), 'multimodel')
                    my_model = obj.study.model.matrixmodel(c, s, b);
                else
                    my_model = obj.study.model;
                end
                nstates = my_model.nstates;
                if nstates>5
                    n1 = ceil(sqrt(nstates));
                    n2 = n1;
                else
                    n1 = 1;
                    n2 = nstates;
                end
                topos = []
                for j=1:nstates
                    subplot(n1, n2, j);
                    topo = my_model.emis_model.posterior.mean_normal{j}.mean;
                    topos = [topos topo];
                    topoplot(topo, ruta);
                    title(['State: ' num2str(j)]);
                end
            end
                
            try
                my_draw(selected(1), selected(2), selected(3));
                if strcmp(class(obj.study.model), 'multimodel')
                    obj.multisubject_selector(f, @my_draw);
                end
            catch err
                obj.my_error(err);
            end
        end
        
        function plot_cov_matrix_action(obj, ~, ~)
            selected = [1, 1, 1];
            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'Covariance Matrix',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            function my_draw(c, s, b)
                if strcmp(class(obj.study.model), 'multimodel')
                    my_model = obj.study.model.matrixmodel(c, s, b);
                else
                    my_model = obj.study.model;
                end
                nstates = my_model.nstates;
                for j=1:nstates
                    subplot(1, nstates, j);
                    cov = my_model.emis_model.parsampl.prec{j};
                    imagesc(cov);
                end
                colormap('jet');
                set(gca,'YDir','normal');
            end

            try
                %David cambia de init(1) a init(3)
                obj.study.model.init(3);
                my_draw(selected(1), selected(2), selected(3));
                if strcmp(class(obj.study.model), 'multimodel')
                    obj.multisubject_selector(f, @my_draw);
                end
            catch err
                obj.my_error(err);
            end
        end
        
        function plot_free_energy(obj, ~, ~)
            figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'Free Energy',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            try
                plot([obj.study.training_results.fedetail.fe]);
                xlabel('Iteration number');
                ylabel('Free Energy');
            catch err
                obj.my_error(err);
            end
        end
        
        function plot_gamma(obj, ~, ~)
            selected = [1, 1, 1];

            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'State Probability',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            function my_draw(c, s, b)
                if strcmp(class(obj.study.model), 'multimodel')
                    results = obj.study.training_results.decodevar.cond(c).subj(s).block(b);
                else
                    results = obj.study.training_results;
                end
                area(results.decodevar.gamma);
                ndata = size(results.decodevar.gamma, 1);
                ylabel('Probability Density');
                xlabel('Time');
                axis([0 ndata 0 1]);
                title('Gamma');
            end
            
            try
                my_draw(selected(1), selected(2), selected(3));
                if obj.study.training_metadata.multisubject
                    obj.multisubject_selector(f, @my_draw);
                    set(f.CurrentAxes, 'Units', 'pixels');
                    set(f.CurrentAxes, 'Position', [f.CurrentAxes.Position(1) f.CurrentAxes.Position(2) f.CurrentAxes.Position(3) f.CurrentAxes.Position(4)-50]);
                end
            catch err
                obj.my_error(err);
            end
        end
        
        function plot_posterior_dur(obj, ~, ~)
            selected = [1, 1, 1];
            is_first = 1;

            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'Duration Posterior',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            function my_draw(c, s, b)
                if ~is_first
                    clf(f.Children(1));
                end
                if strcmp(class(obj.study.model), 'multimodel')
                    my_model = obj.study.model.matrixmodel(c, s, b);
                else
                    my_model = obj.study.model;
                end
                dur = my_model.dur_model.prob('VB',[1:obj.study.opt.dmax]');
                hold on;
                for j=1:size(dur,2)
                    plot(dur(:,j));
                    pal{j}=['State ' num2str(j)];
                end
                xlabel('Duration');
                ylabel('Probability');
                grid;
                legend(pal);
                if obj.study.training_metadata.multisubject && ~is_first
                    set(f.CurrentAxes, 'Units', 'pixels');
                    set(f.CurrentAxes, 'Position', [f.CurrentAxes.Position(1) f.CurrentAxes.Position(2) f.CurrentAxes.Position(3) f.CurrentAxes.Position(4)-50]);
                end
            end
            
            try
                my_draw(selected(1), selected(2), selected(3));
                if obj.study.training_metadata.multisubject
                    obj.multisubject_selector(f, @my_draw);
                    set(f.CurrentAxes, 'Units', 'pixels');
                    set(f.CurrentAxes, 'Position', [f.CurrentAxes.Position(1) f.CurrentAxes.Position(2) f.CurrentAxes.Position(3) f.CurrentAxes.Position(4)-50]);
                    my_panel = f.Children(1);
                end
                is_first = 0;
            catch err
                obj.my_error(err);
            end
        end
        
        function plot_gen_sequence_action(obj, ~, ~)
            selected = [1, 1, 1];

            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'States Sequence',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            try
                my_values = obj.study.generated_metadata.seq;
                b = bar(my_values);
                xlabel('Time');
                ylabel('State');
                title('State Label');
            catch err
                obj.my_error(err);
            end
        end
        
        function plot_gen_occupancy_action(obj, ~, ~)
            selected = [1, 1, 1];
            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'Occupancy',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            try
                if strcmp(class(obj.study.model), 'multimodel')
                    my_model = obj.study.model.matrixmodel(c, s, b);
                else
                    my_model = obj.study.model;
                end
                seq = obj.study.generated_metadata.seq;
                [a d]=hist(seq,1:my_model.nstates);
                bar(d,100.0*a/size(seq,2));
                xlabel('State Label');
                ylabel('% Ocupancy');
                title('State Ocupancy');
            catch err
                obj.my_error(err);
            end
        end
        
        function plot_gen_histogram_action(obj, ~, ~)
            selected = [1, 1, 1];
            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'Duration Histograms',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            try
                if strcmp(class(obj.study.model), 'multimodel')
                    my_model = obj.study.model.matrixmodel(c, s, b);
                else
                    my_model = obj.study.model;
                end
                seq = obj.study.generated_metadata.seq;
                [a d]=util.calculadur(seq);
                for j=1:my_model.nstates
                    subplot(my_model.nstates,1,j);
                    hist(a(d==j));
                    ylabel('Counts');
                    xlabel('State Duration');
                    if j==1
                        title('State Duration Histograms');
                    end
                end
            catch err
                obj.my_error(err);
            end
        end

        function plot_gen_topo_action(obj, ~, ~)
            is_tp = which('topoplot');
            if isempty(is_tp)
                obj.my_warning('It requires Eeglab');
                delete(f);
                return;
            end
            
            ruta = obj.select_eeg_coords();
            if isempty(ruta)
                obj.my_warning('Must choose a coordinates file for your eeg hardware');
                return;
            end
            
            selected = [1, 1, 1];
            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'Topography',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            function my_draw(c, s, b)
                if strcmp(class(obj.study.model), 'multimodel')
                    my_model = obj.study.model.matrixmodel(c, s, b);
                else
                    my_model = obj.study.model;
                end
                nstates = my_model.nstates;
                if nstates>5
                    n1 = ceil(sqrt(nstates));
                    n2 = n1;
                else
                    n1 = 1;
                    n2 = nstates;
                end
                for j=1:nstates
                    subplot(n1, n2, j);
                    topo = my_model.emis_model.parsampl.mean{j}
                    topoplot(topo, ruta);
                    title(['Estado: ' num2str(j)]);
                end
            end
                
            try
                my_draw(selected(1), selected(2), selected(3));
            catch err
                obj.my_error(err);
            end
        end
        
        function plot_gen_cov_matrix_action(obj, ~, ~)
            selected = [1, 1, 1];
            f = figure(...
                'units', 'pixels',...
                'Position', [100 100 600 400],...
                'Name', 'Covariance Matrix',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            try
                if strcmp(class(obj.study.model), 'multimodel')
                    my_model = obj.study.model.matrixmodel(c, s, b);
                else
                    my_model = obj.study.model;
                end
                nstates = my_model.nstates;
                for j=1:nstates
                    subplot(1, nstates, j);
                    cov = my_model.emis_model.parsampl.prec{j};
                    imagesc(cov);
                end
                colormap('jet');
                set(gca,'YDir','normal');
            catch err
                obj.my_error(err);
            end
        end
        
        % dialogs
        function data = get_model(obj)
            data = {};

            model_dialog = figure(...
                'units', 'pixels',...
                'Position', [100 100 350 320],...
                'Toolbar', 'none',...
                'Menubar', 'none',...
                'Name', 'Create Model',...
                'numbertitle', 'off',...
                'resize', 'off');
            % Model HMM o HSMM
            model_selection = uibuttongroup(model_dialog,...
                'Title', 'Model Selection',...
                'units', 'pixels',...
                'Position', [10 260 330 50]);
            model_hsmm = uicontrol(model_selection,...
                'Style', 'radiobutton',...
                'Position', [20 10 110 20],...
                'String', 'HSMM');
            model_hmm = uicontrol(model_selection,...
                'Style', 'radiobutton',...
                'Position', [150 10 110 20],...
                'String', 'HMM');
            % emission dimension
            uicontrol(model_dialog,...
                'Style', 'text',...
                'Position', [20 230 150 20],...
                'String', 'Emission Channels',...
                'HorizontalAlignment', 'left');
            % set the data dimension in the dimension entry
            emis_dim = '';
            if ~isempty(obj.study.data)
                emis_dim = num2str(obj.study.metadata.dimension);
            end
            model_emission_dim_entry = uicontrol(model_dialog,...
                'Style', 'edit',...
                'Position', [145 230 30 20],...
                'String', emis_dim,...
                'HorizontalAlignment', 'left');
            % number of states
            uicontrol(model_dialog,...
                'Style', 'text',...
                'Position', [20 200 150 20],...
                'String', '# States (0: undefined)',...
                'HorizontalAlignment', 'left');
            model_states_entry = uicontrol(model_dialog,...
                'Style', 'edit',...
                'Position', [145 200 30 20],...
                'String', '0',...
                'HorizontalAlignment', 'left');
            % order (for mar model)
            model_order_label = uicontrol(model_dialog,...
                'Style', 'text',...
                'Position', [260 200 30 20],...
                'String', 'Order',...
                'HorizontalAlignment', 'left',...
                'Visible', 'off');
            model_order_entry = uicontrol(model_dialog,...
                'Style', 'edit',...
                'Position', [300 200 30 20],...
                'String', '2',...
                'HorizontalAlignment', 'left',...
                'Visible', 'off');
            
            % models
            model_component = uipanel(model_dialog,...
                'Title', 'Components',...
                'units', 'pixels',...
                'Position', [10 50 330 140]);
            uicontrol(model_component,...
                'Style', 'text',...
                'Position', [10 100 150 20],...
                'String', 'Emission Model',...
                'HorizontalAlignment', 'left');
            emis_opts = obj.classes_ls('emis_model');
            model_emission_combo = uicontrol(model_component,...
                'Style', 'popup',...
                'String', emis_opts,...
                'Position', [160 100 160 20]);
            uicontrol(model_component,...
                'Style', 'text',...
                'Position', [10 70 150 20],...
                'String', 'Transition Model',...
                'HorizontalAlignment', 'left');
            trans_opts = obj.classes_ls('trans_model');
            model_transition_combo = uicontrol(model_component,...
                'Style', 'popup',...
                'String', trans_opts,...
                'Position', [160 70 160 20]);
            uicontrol(model_component,...
                'Style', 'text',...
                'Position', [10 40 150 20],...
                'String', 'Initial Probability Model',...
                'HorizontalAlignment', 'left');
            in_opts = obj.classes_ls('in_model');
            model_initial_combo = uicontrol(model_component,...
                'Style', 'popup',...
                'String', in_opts,...
                'Position', [160 40 160 20]);
            model_duration_label = uicontrol(model_component,...
                'Style', 'text',...
                'Position', [10 10 150 20],...
                'String', 'Duration Model',...
                'HorizontalAlignment', 'left');
            dur_opts = obj.classes_ls('dur_model');
            model_duration_combo = uicontrol(model_component,...
                'Style', 'popup',...
                'String', dur_opts,...
                'Position', [160 10 160 20]);
            % help, create and cancel buttons
            help_button = uicontrol(model_dialog,...
                'Style', 'pushbutton',...
                'String', 'Help',...
                'Position', [10 10 50 20]);
            create_button = uicontrol(model_dialog,...
                'Style', 'pushbutton',...
                'String', 'Create',...
                'ForegroundColor', [0.0 0.8 0.0],...
                'Position', [230 10 50 20]);
            cancel_button = uicontrol(model_dialog,...
                'Style', 'pushbutton',...
                'String', 'Cancel',...
                'ForegroundColor', [0.8 0.0 0.0],...
                'Position', [290 10 50 20]);
            
            function select_emis_combo(~, ~)
                text = emis_opts(model_emission_combo.Value);
                if strcmp(text, 'mar')
                    model_order_label.Visible = 'on';
                    model_order_entry.Visible = 'on';
                else
                    model_order_label.Visible = 'off';
                    model_order_entry.Visible = 'off';
                end
            end
            
            function select_HSMM(~, ~)
                model_duration_label.Visible = 'on';
                model_duration_combo.Visible = 'on';
            end
            
            function select_HMM(~, ~)
                model_duration_label.Visible = 'off';
                model_duration_combo.Visible = 'off';
            end
            
            function help_action(~, ~)
                obj.my_warning('TODO help');
            end
            
            function create_action(~, ~)
                % get info
                selected_model = model_selection.SelectedObject.String;
                dimension = str2num(model_emission_dim_entry.String);
                n_states = str2num(model_states_entry.String);
                
                % validation values are numbers
                if isempty(dimension) || (dimension < 1)
                    obj.my_warning('"Emission Channels" must be integer greater than 0');
                    return;
                end
                if isempty(n_states) ||  (n_states < 0)
                    obj.my_warning('"# Satates" must be integer greater or equal than 0');
                    return;
                end

                % define inner models first
                % emission
                selected_emis = emis_opts(model_emission_combo.Value);
                if strcmp(selected_emis, 'mar')
                    order = str2num(model_order_entry.String);
                    if isempty(order) ||  (order < 0)
                        obj.my_warning('"Order" must be integer greater or equal than 0');
                        return;
                    end
                    selected_emis = eval(['emis_model.' selected_emis{1} '(' num2str(dimension) ',' num2str(n_states)  ',' num2str(order) ')']);
                else
                    selected_emis = eval(['emis_model.' selected_emis{1} '(' num2str(dimension) ',' num2str(n_states) ')']);
                end
                % initial contition
                selected_in = in_opts(model_initial_combo.Value);
                selected_in = eval(['in_model.' selected_in{1} '(' num2str(n_states) ')']);
                % transition
                selected_trans = trans_opts(model_transition_combo.Value);
                selected_trans = eval(['trans_model.' selected_trans{1} '(' num2str(n_states) ')']);
                % duration contition
                if (strcmp(selected_model, 'HSMM'))
                    selected_dur = dur_opts(model_duration_combo.Value);
                    selected_dur = eval(['dur_model.' selected_dur{1} '(1,' num2str(n_states) ')']);
                end
                
                % define model
                if (strcmp(selected_model, 'HSMM'))
                    data.model = hsmm(dimension,...
                        n_states,...
                        selected_emis,...
                        selected_in,...
                        selected_trans,...
                        selected_dur);
                else
                    data.model = hmm(dimension,...
                        n_states,...
                        selected_emis,...
                        selected_in,...
                        selected_trans);
                end
                % close de dialog
                delete(model_dialog);
            end
            
            function cancel_action(~, ~)
                delete(model_dialog);
            end
            
            model_hmm.Callback = @select_HMM;
            model_hsmm.Callback = @select_HSMM;
            model_emission_combo.Callback = @select_emis_combo;
            help_button.Callback = @help_action;
            create_button.Callback = @create_action;
            cancel_button.Callback = @cancel_action;

            select_emis_combo(0, 0);
            uiwait(model_dialog);
        end

        function data = get_distribution(obj, is_prior)
            data = {};
            new_model = struct();
            selected = [1, 1, 1];
            
            function is_equal=compare_structs(elemA,elemB)
                is_equal = 1;
                % if they are not same class, not equal
                if class(elemA) ~= class(elemB)
                    is_equal = 0;
                    return;
                end
                if strcmp(class(elemA), 'struct')
                    % get fields
                    fldsA = sort(fields(elemA));
                    fldsB = sort(fields(elemB));
                    % if they are not same size, not equal
                    if any(any(size(fldsA) ~= size(fldsB)))
                        is_equal = 0;
                        return;
                    end
                    % if subelements are not the same, not equal
                    if ~all(cellfun(@strcmp, fldsA, fldsB))
                        is_equal = 0;
                        return;
                    end
                    % go deeper in subelemets to verify
                    for i=1:size(fldsA,1)
                        for j=1:size(fldsA,2)
                            auxA = getfield(elemA, fldsA{i,j}); 
                            auxB = getfield(elemB, fldsB{i,j});
                            aux_equal = compare_structs(auxA, auxB);
                            if ~aux_equal
                                is_equal = 0;
                                return;
                            end
                        end
                    end
                else
                    if strcmp(class(elemA), 'cell')
                        % if they are not same size, not equal
                        if any(any(size(elemA) ~= size(elemB)))
                            is_equal = 0;
                            return;
                        end
                        % go deeper in subelemets to verify
                        for i=1:size(elemA,1)
                            for j=1:size(elemA,2)
                                aux_equal = compare_structs(elemA{i,j}, elemB{i,j});
                                if ~aux_equal
                                    is_equal = 0;
                                    return;
                                end
                            end
                        end
                    else
                        if strcmp(class(elemA), 'double')
                            is_equal = all(all(elemA == elemB));
                        else
                            is_equal = 0;
                            return;
                        end
                    end
                end
            end
            
            function new_values=non_inf()
                my_model = [];
                if strcmp(class(obj.study.model), 'multimodel')
                    obj.my_warning(['Remember: You cannot re-train a multimodel object,', 'you should create a new model instead']);
                    my_model = obj.study.model.matrixmodel(1,1,1);
                else
                    my_model = obj.study.model;
                end
                model_name = class(my_model);
                % previous values
                previous.emis = my_model.emis_model.prior;
                previous.trans = my_model.trans_model.prior; 
                previous.in = my_model.in_model.prior;
                if strcmp(model_name, 'hsmm')
                    previous.dur = my_model.dur_model.prior;
                end
                % calculate
                my_model.emis_model.priornoinf();
                my_model.trans_model.priornoinf();
                my_model.in_model.priornoinf();
                if strcmp(model_name, 'hsmm')
                    my_model.dur_model.priornoinf();
                end
                % store
                new_values.emis = my_model.emis_model.prior;
                new_values.trans = my_model.trans_model.prior; 
                new_values.in = my_model.in_model.prior;
                if strcmp(model_name, 'hsmm')
                    new_values.dur = my_model.dur_model.prior;
                end
                % set old parameters
                my_model.emis_model.prior = previous.emis;
                my_model.trans_model.prior = previous.trans; 
                my_model.in_model.prior = previous.in;
                if strcmp(model_name, 'hsmm')
                    my_model.dur_model.prior = previous.dur;
                end
            end
            
            function refresh(src, ~, mod, data)
                params = data.ds_params;
                fs = fields(params);
                for i=1:size(fs,1)
                    tb = getfield(params, fs{i});
                    if strcmp(src.String{src.Value}, 'Non Informative')
                        ni = getfield(getfield(noinf, mod), fs{i});
                        tb.UserData = ni;
                        [labels{1:size(tb.Data,1),1:size(tb.Data,2)}] = deal('view');
                    else
                        [labels{1:size(tb.Data,1),1:size(tb.Data,2)}] = deal('edit');
                    end
                    tb.Data = labels;
                    clear('labels');
                end
            end
            
            % non informative data
            if is_prior
                noinf = non_inf();
            end
            
            function cm = copy_model()
                for i=1:size(obj.study.data.cond,2)
                    for j=1:size(obj.study.data.cond(i).subj,2)
                        for k=1:size(obj.study.data.cond(i).subj(j).block,2)
                            aux_m = obj.study.model.matrixmodel(i,j,k);
                            if ~is_prior
                                cm.cond(i).subj(j).block(k).emis = aux_m.emis_model.posterior;
                                if strcmp(class(aux_m), 'hsmm')
                                    cm.cond(i).subj(j).block(k).dur = aux_m.dur_model.posterior;
                                end
                                cm.cond(i).subj(j).block(k).trans = aux_m.trans_model.posterior;
                                cm.cond(i).subj(j).block(k).init = aux_m.in_model.posterior;
                            else
                                cm.cond(i).subj(j).block(k).emis = aux_m.emis_model.prior;
                                if strcmp(class(aux_m), 'hsmm')
                                    cm.cond(i).subj(j).block(k).dur = aux_m.dur_model.prior;
                                end
                                cm.cond(i).subj(j).block(k).trans = aux_m.trans_model.prior;
                                cm.cond(i).subj(j).block(k).init = aux_m.in_model.prior;
                            end
                        end
                    end
                end
            end
            
            % ask if it is multisubject
            if strcmp(class(obj.study.model), 'multimodel')
                multisubject = 1;
                new_model = copy_model();
                my_model = obj.study.model.matrixmodel(selected(1), selected(2), selected(3));
            else
                multisubject = 0;
                if is_prior
                    new_model.emis = obj.study.model.emis_model.prior;
                    if strcmp(class(obj.study.model), 'hsmm')
                        new_model.dur = obj.study.model.dur_model.prior;
                    end
                    new_model.trans = obj.study.model.trans_model.prior;
                    new_model.init = obj.study.model.in_model.prior;
                else
                    new_model.emis = obj.study.model.emis_model.posterior;
                    if strcmp(class(obj.study.model), 'hsmm')
                        new_model.dur = obj.study.model.dur_model.posterior;
                    end
                    new_model.trans = obj.study.model.trans_model.posterior;
                    new_model.init = obj.study.model.in_model.posterior;
                end
                my_model = obj.study.model;
            end
            % global model name
            gm_name = class(my_model);
            
            dist_dialog = figure(...
                'units', 'pixels',...
                'Position', [100 50 500 400],...
                'Toolbar', 'none',...
                'Menubar', 'none',...
                'Name', 'Prior Distributions',...
                'numbertitle', 'off',...
                'resize', 'off',...
                'UserData', my_model);
            if ~is_prior
                dist_dialog.Name = 'Posterior Distributions';
            end

            x = 0;
            y = 40;
            % THIS ONE DEPENDS ON MODEL (HMM or HSMM)
            % duration panel
            duration_panel = {};
            if(strcmp(gm_name, 'hsmm'))
                duration_panel = obj.get_emission_duration_panel(dist_dialog, 'Duration', my_model.dur_model, is_prior);
                aux = duration_panel.ds_panel.Position;
                duration_panel.ds_panel.Position = [10 y aux(3) aux(4)];
                if is_prior
                    duration_panel.ds_combo.Callback = {@refresh, 'dur', duration_panel};
                    if ~compare_structs(noinf.dur, my_model.dur_model.prior)
                        duration_panel.ds_combo.Value = 2;
                        refresh(duration_panel.ds_combo, 0, 'dur', duration_panel);
                    end
                end
                x = max(x, aux(3));
                y = y + 10 + aux(4);
            end
            % initial condition panel
            init_panel = obj.get_transition_initial_panel(dist_dialog, 'Initial Condition', my_model.in_model, is_prior);
            aux = init_panel.ds_panel.Position;
            init_panel.ds_panel.Position = [10 y aux(3) aux(4)];
            if is_prior
                init_panel.ds_combo.Callback = {@refresh, 'in', init_panel};
                if ~compare_structs(noinf.in, my_model.in_model.prior)
                    init_panel.ds_combo.Value = 2;
                    refresh(init_panel.ds_combo, 0, 'dur', init_panel);
                end
            end
            x = max(x, aux(3));
            y = y + 10 + aux(4);
            % transition panel
            transition_panel = obj.get_transition_initial_panel(dist_dialog, 'Transition', my_model.trans_model, is_prior);
            aux = transition_panel.ds_panel.Position;
            transition_panel.ds_panel.Position = [10 y aux(3) aux(4)];
            if is_prior
                transition_panel.ds_combo.Callback = {@refresh, 'trans', transition_panel};
                if ~compare_structs(noinf.trans, my_model.trans_model.prior)
                    transition_panel.ds_combo.Value = 2;
                    refresh(transition_panel.ds_combo, 0, 'dur', transition_panel);
                end
            end
            x = max(x, aux(3));
            y = y + 10 + aux(4);
            % emission panel
            emission_panel = obj.get_emission_duration_panel(dist_dialog, 'Emission', my_model.emis_model, is_prior);
            aux = emission_panel.ds_panel.Position;
            emission_panel.ds_panel.Position = [10 y aux(3) aux(4)];
            if is_prior
                emission_panel.ds_combo.Callback = {@refresh, 'emis', emission_panel};
                if ~compare_structs(noinf.emis, my_model.emis_model.prior)
                    emission_panel.ds_combo.Value = 2;
                    refresh(emission_panel.ds_combo, 0, 'dur', emission_panel);
                end
            end
            x = max(x, aux(3))+20;
            y = y + 10 + aux(4);
            
            % recalculate the size of the dialog
            dist_dialog.Position = [100 50 x y];
            
            if(strcmp(gm_name, 'hsmm'))
                duration_panel.ds_panel.Position(3) = x-20;
            end
            init_panel.ds_panel.Position(3) = x-20;
            transition_panel.ds_panel.Position(3) = x-20;
            emission_panel.ds_panel.Position(3) = x-20;

            % help, create and cancel buttons
            help_button = uicontrol(dist_dialog,...
                'Style', 'pushbutton',...
                'String', 'Help',...
                'Position', [10 10 50 20]);
            
            if is_prior
                save_button = uicontrol(dist_dialog,...
                    'Style', 'pushbutton',...
                    'String', 'Save',...
                    'ForegroundColor', [0.0 0.8 0.0],...
                    'Position', [x-120 10 50 20]);
                cancel_button = uicontrol(dist_dialog,...
                    'Style', 'pushbutton',...
                    'String', 'Cancel',...
                    'ForegroundColor', [0.8 0.0 0.0],...
                    'Position', [x-60 10 50 20]);
            else
                ok_button = uicontrol(dist_dialog,...
                    'Style', 'pushbutton',...
                    'String', 'OK',...
                    'ForegroundColor', [0.0 0.8 0.0],...
                    'Position', [x-60 10 50 20]);
            end
            
            function save_model()
                model_name = '';
                if multisubject
                    model_name = class(obj.study.model.matrixmodel(selected(1), selected(2), selected(3)));
                else
                    model_name = class(obj.study.model);
                end
                % emission data
                emis_params = {};
                f1 = fields(emission_panel.ds_params);
                for i=1:size(f1,1)
                    aux1 = getfield(emission_panel.ds_params, f1{i});
                    aux2 = aux1.UserData;
                    emis_params = setfield(emis_params, f1{i}, aux2);
                end
                % duration data
                if strcmp(model_name, 'hsmm')
                    dur_params = {};
                    f1 = fields(duration_panel.ds_params);
                    for i=1:size(f1,1)
                        aux1 = getfield(duration_panel.ds_params, f1{i});
                        aux2 = aux1.UserData;
                        dur_params = setfield(dur_params, f1{i}, aux2);
                    end
                end
                % transition data
                trans_params = {};
                f1 = fields(transition_panel.ds_params);
                for i=1:size(f1,1)
                    aux1 = getfield(transition_panel.ds_params, f1{1});
                    aux2 = aux1.UserData;
                    trans_params = setfield(trans_params, f1{1}, aux2);
                end
                % initial data
                init_params = {};
                f1 = fields(init_panel.ds_params);
                for i=1:size(f1,1)
                    aux1 = getfield(init_panel.ds_params, f1{1});
                    aux2 = aux1.UserData;
                    init_params = setfield(init_params, f1{1}, aux2);
                end
                % put data in the models
                if multisubject
                    new_model.cond(selected(1)).subj(selected(2)).block(selected(3)).emis = emis_params;
                    if strcmp(model_name, 'hsmm')
                        new_model.cond(selected(1)).subj(selected(2)).block(selected(3)).dur = dur_params;
                    end
                    new_model.cond(selected(1)).subj(selected(2)).block(selected(3)).trans = trans_params;
                    new_model.cond(selected(1)).subj(selected(2)).block(selected(3)).init = init_params;
                else
                    new_model.emis = emis_params;
                    if strcmp(model_name, 'hsmm')
                        new_model.dur = dur_params;
                    end
                    new_model.trans = trans_params;
                    new_model.init = init_params;
                end
            end
            
            function my_callback(i_cond, i_subj, i_block)
                if is_prior
                    save_model();
                end
                selected = [i_cond i_subj i_block];
                % emission
                fs = fields(emission_panel.ds_params);
                for i=1:size(fs, 1)
                    aux_val = getfield(new_model.cond(selected(1)).subj(selected(2)).block(selected(3)).emis, fs{1});
                    ud = getfield(emission_panel.ds_params, fs{1});
                    ud.UserData = aux_val;
                    emission_panel.ds_params = setfield(emission_panel.ds_params, fs{1}, ud);
                end
                % duration
                if strcmp(class(obj.study.model.matrixmodel(1,1,1)), 'hsmm')
                    fs = fields(duration_panel.ds_params);
                    for i=1:size(fs, 1)
                        aux_val = getfield(new_model.cond(selected(1)).subj(selected(2)).block(selected(3)).dur, fs{1});
                        ud = getfield(duration_panel.ds_params, fs{1});
                        ud.UserData = aux_val;
                        duration_panel.ds_params = setfield(duration_panel.ds_params, fs{1}, ud);
                    end
                end
                % initial condition
                fs = fields(init_panel.ds_params);
                for i=1:size(fs, 1)
                    aux_val = getfield(new_model.cond(selected(1)).subj(selected(2)).block(selected(3)).init, fs{1});
                    ud = getfield(init_panel.ds_params, fs{1});
                    ud.UserData = aux_val;
                    init_panel.ds_params = setfield(init_panel.ds_params, fs{1}, ud);
                end
                % transition
                fs = fields(transition_panel.ds_params);
                for i=1:size(fs, 1)
                    aux_val = getfield(new_model.cond(selected(1)).subj(selected(2)).block(selected(3)).trans, fs{1});
                    ud = getfield(transition_panel.ds_params, fs{1});
                    ud.UserData = aux_val;
                    transition_panel.ds_params = setfield(transition_panel.ds_params, fs{1}, ud);
                end
            end
            
            if multisubject && ~is_prior
                obj.multisubject_selector(dist_dialog, @my_callback);
            end
            
            function help_action(~, ~)
                obj.my_warning('TODO help');
            end
            
            function save_action(~, ~)
                save_model();
                data.model = new_model;
                data.multisubject = multisubject;
                delete(dist_dialog);
            end
            
            function cancel_action(~, ~)
                delete(dist_dialog);
            end
            
            help_button.Callback = @help_action;
            
            if is_prior
                save_button.Callback = @save_action;
                cancel_button.Callback = @cancel_action;
            else
                ok_button.Callback = @cancel_action;
            end

            data = {};
            
            uiwait(dist_dialog);
        end
        
        function data = get_estimation(obj)
            data = {};
            multisubject = 0;
            
            est_dialog = figure(...
                'units', 'pixels',...
                'Position', [100 100 240 330],...
                'Toolbar', 'none',...
                'Menubar', 'none',...
                'Name', 'Train Model',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            y = 0;
            if obj.study.metadata.multisubject
                multisubject = 1;
                p = uipanel(est_dialog,...
                    'units', 'pixels',...
                    'Title', '',...
                    'FontWeight', 'bold',...
                    'Position', [30 40 180 100]);
                % block
                uicontrol(p, 'Style', 'text', 'String', 'block', 'HorizontalAlignment', 'right', 'Position', [5 5 35 20]);
                check.emis(3) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [48 8 30 20]);
                check.trans(3) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [83 8 30 20]);
                check.in(3) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [118 8 30 20]);
                check.dur(3) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [153 8 30 20]);
                % subject
                uicontrol(p, 'Style', 'text', 'String', 'subj', 'HorizontalAlignment', 'right', 'Position', [5 30 35 20]);
                check.emis(2) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [48 33 30 20]);
                check.trans(2) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [83 33 30 20]);
                check.in(2) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [118 33 30 20]);
                check.dur(2) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [153 33 30 20]);
                % condition
                uicontrol(p, 'Style', 'text', 'String', 'cond', 'HorizontalAlignment', 'right', 'Position', [5 55 35 20]);
                check.emis(1) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [48 58 30 20]);
                check.trans(1) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [83 58 30 20]);
                check.in(1) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [118 58 30 20]);
                check.dur(1) = uicontrol(p, 'Style', 'checkbox', 'String', '', 'HorizontalAlignment', 'center', 'Position', [153 58 30 20]);
                % header
                uicontrol(p, 'Style', 'text', 'String', 'emis', 'HorizontalAlignment', 'center', 'Position', [40 75 30 20]);
                uicontrol(p, 'Style', 'text', 'String', 'trans', 'HorizontalAlignment', 'center', 'Position', [75 75 30 20]);
                uicontrol(p, 'Style', 'text', 'String', 'in', 'HorizontalAlignment', 'center', 'Position', [110 75 30 20]);
                uicontrol(p, 'Style', 'text', 'String', 'dur', 'HorizontalAlignment', 'center', 'Position', [145 75 30 20]);
                y = y+100;
            end
            
            % learning algorithm
            uicontrol(est_dialog,...
                'Style', 'text',...
                'String', 'Learning Algorithm:',...
                'HorizontalAlignment', 'left',...
                'Position', [10 y+295 100 20]);
            alg_idx = find(strcmp({'VB', 'EM'}, obj.study.opt.train));
            alg_opts = {'Variational Bayes', 'Expectation Maximization'};
            algorithm_combo = uicontrol(est_dialog,...
                'Style', 'popup',...
                'String', alg_opts,...
                'Position', [110 y+300 120 20],...
                'Value', alg_idx);
            % Initialization
            is_kmeans = 0;
            init_panel = uipanel(est_dialog,...
                'Title', 'Initialization',...
                'units', 'pixels',...
                'Position', [10 y+225 220 70]);
            uicontrol(init_panel,...
                'Style', 'text',...
                'String', 'Algorithm:',...
                'HorizontalAlignment', 'left',...
                'Position', [5 28 100 20]);
            init_idx = find(strcmp({'kmeans', 'random'}, obj.study.opt.train));
            init_opts = {'K-means', 'Random'};
            init_combo = uicontrol(init_panel,...
                'Style', 'popup',...
                'String', init_opts,...
                'Position', [100 30 110 20]);
            % max iter kmeans
            km_lab = uicontrol(init_panel,...
                'Style', 'text',...
                'String', 'Iterations:',...
                'HorizontalAlignment', 'left',...
                'Position', [220 28 50 20]);
            km_entry = uicontrol(init_panel,...
                'Style', 'edit',...
                'String', '',...
                'Position', [280 30 40 20],...
                'String', num2str(obj.study.opt.maxiter));
            % repetitions init
            repetitions_label = uicontrol(init_panel,...
                'Style', 'text',...
                'String', 'Repetitions:',...
                'HorizontalAlignment', 'left',...
                'Position', [5 2 60 20]);
            repetitions_entry = uicontrol(init_panel,...
                'Style', 'edit',...
                'String', '',...
                'Position', [100 5 110 20],...
                'HorizontalAlignment', 'right',...
                'String', num2str(obj.study.opt.nrep));
            % tolerance
            uicontrol(est_dialog,...
                'Style', 'text',...
                'String', 'Tolerance:',...
                'HorizontalAlignment', 'left',...
                'Position', [10 y+200 70 20]);
            tolerance_entry = uicontrol(est_dialog,...
                'Style', 'edit',...
                'Position', [110 y+200 120 20],...
                'HorizontalAlignment', 'right',...
                'String', num2str(obj.study.opt.tol));
            % max iter global
            uicontrol(est_dialog,...
                'Style', 'text',...
                'String', 'Iterations:',...
                'HorizontalAlignment', 'left',...
                'Position', [10 y+170 70 20]);
            iterations_entry = uicontrol(est_dialog,...
                'Style', 'edit',...
                'Position', [110 y+170 120 20],...
                'HorizontalAlignment', 'right',...
                'String', num2str(obj.study.opt.maxcyc));
            % repetitions global
            uicontrol(est_dialog,...
                'Style', 'text',...
                'String', 'Repetitions:',...
                'HorizontalAlignment', 'left',...
                'Position', [10 y+140 70 20]);
            repetitions_global_entry = uicontrol(est_dialog,...
                'Style', 'edit',...
                'Position', [110 y+140 120 20],...
                'HorizontalAlignment', 'right',...
                'String', 1);
            % duration
            uicontrol(est_dialog,...
                'Style', 'text',...
                'String', 'Duration:',...
                'HorizontalAlignment', 'left',...
                'Position', [10 y+110 70 20]);
            dmax_entry = uicontrol(est_dialog,...
                'Style', 'edit',...
                'Position', [110 y+110 120 20],...
                'HorizontalAlignment', 'right',...
                'String', num2str(obj.study.opt.dmax));
            % # states
            states_label = uicontrol(est_dialog,...
                'Style', 'text',...
                'String', '# States:',...
                'HorizontalAlignment', 'left',...
                'Position', [10 y+80 50 20]);
            states_min_label = uicontrol(est_dialog,...
                'Style', 'text',...
                'String', 'min: ',...
                'HorizontalAlignment', 'right',...
                'Position', [80 y+80 30 20]);
            min_entry = uicontrol(est_dialog,...
                'Style', 'edit',...
                'Position', [110 y+80 40 20],...
                'HorizontalAlignment', 'right',...
                'String', num2str(obj.study.opt.minstates));
            states_max_label = uicontrol(est_dialog,...
                'Style', 'text',...
                'String', 'max: ',...
                'HorizontalAlignment', 'right',...
                'Position', [160 y+80 30 20]);
            max_entry = uicontrol(est_dialog,...
                'Style', 'edit',...
                'Position', [190 y+80 40 20],...
                'HorizontalAlignment', 'right',...
                'String', num2str(obj.study.opt.maxstates));
            if obj.study.model.nstates > 0
               states_label.Visible = 'off';
               states_min_label.Visible = 'off';
               min_entry.Visible = 'off';
               states_max_label.Visible = 'off';
               max_entry.Visible = 'off';
            end
            verbose_check = uicontrol(est_dialog,...
                'Style', 'checkbox',...
                'String', 'Verbose',...
                'Position', [10 y+50 100 20],...
                'Value', strcmp(obj.study.opt.verbose, 'yes'));
            
            % calculate the new high
            est_dialog.Position(4) = est_dialog.Position(4) + y;
            
            % help, create and cancel buttons
            help_button = uicontrol(est_dialog,...
                'Style', 'pushbutton',...
                'String', 'Help',...
                'Position', [10 10 50 20]);
            run_button = uicontrol(est_dialog,...
                'Style', 'pushbutton',...
                'String', 'Run',...
                'ForegroundColor', [0.0 0.8 0.0],...
                'Position', [120 10 50 20]);
            cancel_button = uicontrol(est_dialog,...
                'Style', 'pushbutton',...
                'String', 'Cancel',...
                'ForegroundColor', [0.8 0.0 0.0],...
                'Position', [180 10 50 20]);
            
            function kmeans_options(src, ev)
                if strcmp(init_opts(init_combo.Value), 'K-means')
                    if ~is_kmeans
                        km_lab.Visible = 'on';
                        km_entry.Visible = 'on';
                        run_button.Position(1) = run_button.Position(1) + 110;
                        cancel_button.Position(1) = cancel_button.Position(1) + 110;
                        est_dialog.Position(3) = est_dialog.Position(3) + 110;
                        init_panel.Position(3) = init_panel.Position(3) + 110;
                        repetitions_label.Enable = 'on';
                        repetitions_entry.Enable = 'on';
                        is_kmeans = 1;
                    end
                else
                    if is_kmeans
                        km_lab.Visible = 'off';
                        km_entry.Visible = 'off';
                        run_button.Position(1) = run_button.Position(1) - 110;
                        cancel_button.Position(1) = cancel_button.Position(1) - 110;
                        est_dialog.Position(3) = est_dialog.Position(3) - 110;
                        init_panel.Position(3) = init_panel.Position(3) - 110;
                        repetitions_label.Enable = 'off';
                        repetitions_entry.Enable = 'off';
                        is_kmeans = 0;
                    end
                end
            end
            
            function help_action(src, ev)
                obj.my_warning('TODO help');
            end
            
            function run_action(src, ev)
                % get info
                selected_alg = alg_opts(algorithm_combo.Value);
                selected_init = init_opts(init_combo.Value);
                max_iter = str2num(km_entry.String);
                if isempty(max_iter) || (max_iter < 1)
                    obj.my_warning('"Iterations" (k-means) must be integer greater than 0');
                    return;
                end
                n_repetitions = str2num(repetitions_entry.String);
                if isempty(n_repetitions) || (n_repetitions < 1)
                    obj.my_warning('"Repetitions" (initialization algorithm) must be integer greater than 0');
                    return;
                end
                tolerance = str2num(tolerance_entry.String);
                if isempty(tolerance) || (tolerance <= 0)
                    obj.my_warning('"Tolerance" must be greater than 0');
                    return;
                end
                maxitersim = str2num(iterations_entry.String);
                if isempty(maxitersim) || (maxitersim < 1)
                    obj.my_warning('"Iterations" (global algorithm) must be integer greater than 0');
                    return;
                end
                max_cycle = str2num(repetitions_global_entry.String);
                if isempty(max_cycle) || (max_cycle < 1)
                    obj.my_warning('"Repetitions" (global algorithm) must be integer greater than 0');
                    return;
                end
                dmax = str2num(dmax_entry.String);
                if isempty(dmax) || (dmax < 1)
                    obj.my_warning('"Duration" must be integer greater than 0');
                    return;
                end
                min_states = str2num(min_entry.String);
                if isempty(min_states) || (min_states < 1)
                    obj.my_warning('"States Min" must be integer greater than 0');
                    return;
                end
                max_states = str2num(max_entry.String);
                if isempty(max_states) || (max_states < 1) || (max_states < min_states)
                    obj.my_warning('"States Max" must be integer greater than 0 and greater or equal than "States Min"');
                    return;
                end
                verbose = verbose_check.Value;
                
                % validate values are numbers
                
                % store info
                data.learning_algorithm = selected_alg{1};
                data.initialization = selected_init{1};
                data.maxiter = max_iter;
                data.n_repetitions = n_repetitions;
                data.tolerance= tolerance;
                data.max_iter_sim = max_cycle;
                data.max_cycles = maxitersim;
                data.dmax = dmax;
                data.n_states = {min_states, max_states};
                data.verbose = verbose;
                
                if multisubject
                    data.shareemis = [check.emis(1).Value, check.emis(2).Value, check.emis(3).Value];
                    data.sharein = [check.in(1).Value, check.in(2).Value, check.in(3).Value];
                    data.sharetrans = [check.trans(1).Value, check.trans(2).Value, check.trans(3).Value];
                    data.sharedur = [check.dur(1).Value, check.dur(2).Value, check.dur(3).Value];
                end
                
                data.multisubject = multisubject;
                % close the dialog
                delete(est_dialog);
            end
            
            function cancel_action(src, ev)
                delete(est_dialog);
            end
            
            kmeans_options({},{});
            init_combo.Callback = @kmeans_options;
            help_button.Callback = @help_action;
            run_button.Callback = @run_action;
            cancel_button.Callback = @cancel_action;

            uiwait(est_dialog);
        end
        
        function data = get_generation(obj)
            model_name = class(obj.study.model);
            editable_params = 1;
            mean_based_params = {};
            n_states = obj.study.model.nstates;
            
            function new_values=get_params(option)
                % previous values
                previous.emis = obj.study.model.emis_model.parsampl;
                previous.trans = obj.study.model.trans_model.parsampl; 
                previous.in = obj.study.model.in_model.parsampl;
                if strcmp(model_name, 'hsmm')
                    previous.dur = obj.study.model.dur_model.parsampl;
                end
                % calculate
                obj.study.model.init(option);
                % store
                new_values.emis = obj.study.model.emis_model.parsampl;
                new_values.trans = obj.study.model.trans_model.parsampl; 
                new_values.in = obj.study.model.in_model.parsampl;
                if strcmp(model_name, 'hsmm')
                    new_values.dur = obj.study.model.dur_model.parsampl;
                end
                % get distribution model names
                new_values.models.emis = class(obj.study.model.emis_model);
                new_values.models.emis = new_values.models.emis(strfind(new_values.models.emis, '.')+1:size(new_values.models.emis, 2));
                new_values.models.trans = class(obj.study.model.trans_model);
                new_values.models.trans = new_values.models.trans(strfind(new_values.models.trans, '.')+1:size(new_values.models.trans, 2));
                new_values.models.in = class(obj.study.model.in_model);
                new_values.models.in = new_values.models.in(strfind(new_values.models.in, '.')+1:size(new_values.models.in, 2));
                if strcmp(model_name, 'hsmm')
                    new_values.models.dur = class(obj.study.model.dur_model);
                    new_values.models.dur = new_values.models.dur(strfind(new_values.models.dur, '.')+1:size(new_values.models.dur, 2));
                end
                % set old parameters
                obj.study.model.emis_model.parsampl = previous.emis;
                obj.study.model.trans_model.parsampl = previous.trans; 
                obj.study.model.in_model.parsampl = previous.in;
                if strcmp(model_name, 'hsmm')
                    obj.study.model.dur_model.parsampl = previous.dur;
                end
            end
            
            user_defined_params = get_params(4);

            function modify_cell(src, ev)
                % nobody generate this event so skip
                if isempty(ev.Indices)
                    return
                end
                
                col = ev.Indices(2); % n_states
                row = ev.Indices(1); % parameter name
                
                text = 'weiv';
                if editable_params
                    text = 'tide';
                end
                [aux{1:size(src.Data,1),1:size(src.Data,2)}] = deal(text);
                src.Data = aux;
                text = 'view';
                if editable_params
                    text = 'edit';
                end
                [aux{1:size(src.Data,1),1:size(src.Data,2)}] = deal(text);
                src.Data = aux;

                row_name = src.RowName{row};
                the_model = src.UserData;
                sub_dialog = figure(...
                    'units', 'pixels',...
                    'Position', [100 100 182 166],...
                    'Toolbar', 'none',...
                    'Menubar', 'none',...
                    'Name', '',...
                    'numbertitle', 'off',...
                    'resize', 'off');
               
                f1 = getfield(the_model, row_name);
                f1 = f1{col};
                cw = num2cell(40*ones(1, size(f1, 2)));
                ed = false(1, size(f1, 2));
                if editable_params
                    ed = true(1, size(f1, 2));
                end

                my_table = uitable(sub_dialog,...
                    'ColumnEditable', ed,...
                    'Data', f1,...
                    'ColumnWidth', cw,...
                    'Position', [10 40 162 86]);

                the_label = uicontrol(sub_dialog,...
                    'Position', [10 136 162 20],...
                    'Style', 'text',...
                    'HorizontalAlignment', 'left',...
                    'String', ['State N°' num2str(col) ' / parameter: ' row_name]);

                function load_action(~, ~)
                    % retrieve var from global workspace
                    vars = evalin('base', 'who()');
                    % ask the user to select one
                    [indx,tf] = listdlg('PromptString', 'Select Data Variable:',...
                        'ListSize', [250 150],...
                        'SelectionMode', 'single',...
                        'ListString', vars);
                    if ~tf
                        % cancel the loading by close or press cancel
                        return;
                    end
                    % get the value
                    value = evalin('base', vars{indx});
                    my_table.Data = value;
                end
                    
                function ok_action(~, ~)
                    isn = find(isnan(my_table.Data));
                    if ~isempty(isn) && size(isn, 1) > 0
                        obj.my_warning('Matrix (or vector) cannot contains NaN values');
                        return;
                    end
                    new_data = getfield(the_model, row_name);
                    new_data{col} = my_table.Data;
                    the_model = setfield(the_model, row_name, new_data);
                    src.UserData = the_model;
                    delete(sub_dialog);
                end
                
                function cancel_action(~, ~)
                    delete(sub_dialog);
                end

                if editable_params
                    load_button = uicontrol(sub_dialog,...
                        'Position', [122 136 50 20],...
                        'Style', 'pushbutton',...
                        'String', 'load');

                    load_button.Callback = @load_action;
                    
                    the_label.Position(2) = the_label.Position(2)+30;
                    sub_dialog.Position(4) = sub_dialog.Position(4)+30;
                    
                    ok_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'OK',...
                        'ForegroundColor', [0.0 0.8 0.0],...
                        'Position', [62 10 50 20]);
                    cancel_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'Cancel',...
                        'ForegroundColor', [0.8 0.0 0.0],...
                        'Position', [122 10 50 20]);
                    ok_button.Callback = @ok_action;
                    cancel_button.Callback = @cancel_action;
                else
                    ok_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'OK',...
                        'ForegroundColor', [0.0 0.8 0.0],...
                        'Position', [122 10 50 20]);
                    ok_button.Callback = @cancel_action;
                end
            end
                            
            function modify_trans_in(src, ev)
                the_model = getfield(user_defined_params, src.UserData);
                name = 'Initial Probabilitiy';
                if strcmp(src.UserData, 'trans')
                    name = 'Transitions Matrix';
                end
                    
                sub_dialog = figure(...
                    'units', 'pixels',...
                    'Position', [100 100 182 166],...
                    'Toolbar', 'none',...
                    'Menubar', 'none',...
                    'Name', name,...
                    'numbertitle', 'off',...
                    'resize', 'off');
                
                the_label = uicontrol(sub_dialog,...
                    'Position', [10 136 162 20],...
                    'Style', 'text',...
                    'HorizontalAlignment', 'left',...
                    'String', name);

                cw = num2cell(40*ones(1, size(the_model, 2)));
                ed = true(1, size(the_model, 2));
                if ~editable_params
                    ed = false(1, size(the_model, 2));
                end

                my_table = uitable(sub_dialog,...
                    'ColumnEditable', ed,...
                    'Data', the_model,...
                    'ColumnWidth', cw,...
                    'Position', [10 40 162 86]);

                function ok_action(~, ~)
                    isn = find(isnan(my_table.Data));
                    if ~isempty(isn) && size(isn, 1) > 0
                        obj.my_warning('Matrix (or vector) cannot contains NaN values');
                        return;
                    end
                    user_defined_params = setfield(user_defined_params, src.UserData, my_table.Data);
                    delete(sub_dialog);
                end
                
                function cancel_action(~, ~)
                    delete(sub_dialog);
                end

                function load_action(~, ~)
                    % retrieve var from global workspace
                    vars = evalin('base', 'who()');
                    % ask the user to select one
                    [indx,tf] = listdlg('PromptString', 'Select Data Variable:',...
                        'SelectionMode', 'single',...
                        'ListString', vars);
                    if ~tf
                        % cancel the loading by close or press cancel
                        return;
                    end
                    % get the value
                    value = evalin('base', vars{indx});
                    my_table.Data = value;
                end
                    
                if editable_params
                        load_button = uicontrol(sub_dialog,...
                        'Position', [122 136 50 20],...
                        'Style', 'pushbutton',...
                        'String', 'load');

                    load_button.Callback = @load_action;
                    
                    the_label.Position(2) = the_label.Position(2)+30;
                    sub_dialog.Position(4) = sub_dialog.Position(4)+30;

                    ok_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'OK',...
                        'ForegroundColor', [0.0 0.8 0.0],...
                        'Position', [62 10 50 20]);
                    cancel_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'Cancel',...
                        'ForegroundColor', [0.8 0.0 0.0],...
                        'Position', [122 10 50 20]);
                    ok_button.Callback = @ok_action;
                    cancel_button.Callback = @cancel_action;
                else
                    ok_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'OK',...
                        'ForegroundColor', [0.0 0.8 0.0],...
                        'Position', [122 10 50 20]);
                    ok_button.Callback = @cancel_action;
                end
            end
                            
            data = {};

            gen_dialog = figure(...
                'units', 'pixels',...
                'Position', [100 100 350 480],...
                'Toolbar', 'none',...
                'Menubar', 'none',...
                'Name', 'Generate Data',...
                'numbertitle', 'off',...
                'resize', 'off');

            y = 40;
            if strcmp(model_name, 'hsmm')
                % duration model data 
                f_dur = fields(obj.study.model.dur_model.parsampl);
                [dur_data{1:size(f_dur,1),1:n_states}] = deal('edit');
                [cw{1:n_states}] = deal(30);
                duration_panel = uipanel(gen_dialog,...
                    'Title', 'Duration Model',...
                    'units', 'pixels',...
                    'Position', [10 y 330 140]);
                duration_model_label = uicontrol(duration_panel,...
                    'Position', [10 100 300 20],...
                    'Style', 'text',...
                    'HorizontalAlignment', 'left',...
                    'String', user_defined_params.models.dur);
                duration_table = uitable(duration_panel,...
                    'RowName', f_dur,...
                    'Data', dur_data,...
                    'ColumnWidth', cw,...
                    'Position', [10 10 310 90],...
                    'CellSelectionCallback', @modify_cell,...
                    'UserData', user_defined_params.dur);
            y = y + 150;
            end
            
            % emission model data 
            f_emis = fields(obj.study.model.emis_model.parsampl);
            [emis_data{1:size(f_emis,1),1:n_states}] = deal('edit');
            [cw{1:n_states}] = deal(30);
            emission_panel = uipanel(gen_dialog,...
                'Title', 'Emission Model',...
                'units', 'pixels',...
                'Position', [10 y 330 140]);
            emission_model_label = uicontrol(emission_panel,...
                'Position', [10 100 300 20],...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'String', user_defined_params.models.emis);
            emis_table = uitable(emission_panel,...
                'RowName', {'mean', 'prec'},...
                'Data', emis_data,...
                'ColumnWidth', cw,...
                'Position', [10 10 310 90],...
                'CellSelectionCallback', @modify_cell,...
                'UserData', user_defined_params.emis);

            y = y + 150;
            % transition matrix
            transition_label = uicontrol(gen_dialog,...
                'Position', [20 y 60 40],...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'String', 'Transition Matrix');
            transition_button = uicontrol(gen_dialog,...
                'Style', 'pushbutton',...
                'String', 'edit',...
                'Position', [90 y+20 60 20],...
                'Callback', @modify_trans_in,...
                'UserData', 'trans');
            
            % initial conditions
            initial_label = uicontrol(gen_dialog,...
                'Position', [170 y 60 40],...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'String', 'Initial Probability');
            initial_button = uicontrol(gen_dialog,...
                'Style', 'pushbutton',...
                'String', 'edit',...
                'Position', [240 y+20 60 20],...
                'Callback', @modify_trans_in,...
                'UserData', 'in');
            
            y = y + 50;
            % parameters 
            generation_based_on = uibuttongroup(gen_dialog,...
                'Title', 'Parameters',...
                'units', 'pixels',...
                'Position', [10 y 300 50]);
            user_defined = uicontrol(generation_based_on,...
                'Style', 'radiobutton',...
                'Position', [20 10 110 20],...
                'String', 'User Defined');
            mean_based = uicontrol(generation_based_on,...
                'Style', 'radiobutton',...
                'Position', [150 10 130 20],...
                'String', 'Posterior Parameters');
            if ~obj.study.model.emis_model.posteriorfull()
               mean_based.Enable = 'off'; 
            end

            y = y + 60;
            % data amount
            uicontrol(gen_dialog,...
                'Position', [10 y-3 60 20],...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'String', 'Data Length');
            data_amount_entry = uicontrol(gen_dialog,...
                'Position', [80 y 30 20],...
                'Style', 'edit',...
                'HorizontalAlignment', 'left');
            
            y = y + 30;
            gen_dialog.Position(4) = y;
            
            % help, create and cancel buttons
            help_button = uicontrol(gen_dialog,...
                'Style', 'pushbutton',...
                'String', 'Help',...
                'Position', [10 10 50 20]);
            generate_button = uicontrol(gen_dialog,...
                'Style', 'pushbutton',...
                'String', 'Generate',...
                'ForegroundColor', [0.0 0.8 0.0],...
                'Position', [230 10 50 20]);
            cancel_button = uicontrol(gen_dialog,...
                'Style', 'pushbutton',...
                'String', 'Cancel',...
                'ForegroundColor', [0.8 0.0 0.0],...
                'Position', [290 10 50 20]);

            function mean_based_action(src, ~)
                % ask if model created and if it was trained
                if isempty(mean_based_params)
                    mean_based_params = get_params(1);
                end
                user_defined_params = mean_based_params;
                editable_params = 0;
                [aux1{1:size(emis_table.Data,1),1:size(emis_table.Data,2)}] = deal('view');
                emis_table.Data = aux1;
                [aux2{1:size(duration_table.Data,1),1:size(duration_table.Data,2)}] = deal('view');
                duration_table.Data = aux2;
                transition_button.String = 'view';
                initial_button.String = 'view';
            end
            
            function user_defined_action(~, ~)
                editable_params = 1;
                [aux1{1:size(emis_table.Data,1),1:size(emis_table.Data,2)}] = deal('edit');
                emis_table.Data = aux1;
                [aux2{1:size(duration_table.Data,1),1:size(duration_table.Data,2)}] = deal('edit');
                duration_table.Data = aux2;
                transition_button.String = 'edit';
                initial_button.String = 'edit';
            end
            
            function help_action(~, ~)
                obj.my_warning('TODO help');
            end
            
            function generate_action(~, ~)
                % get info
                selected_model = generation_based_on.SelectedObject.String;
                n_data = str2num(data_amount_entry.String);
                % validation values are numbers
                if isempty(n_data) || (n_data < 1)
                    obj.my_warning('"Data Size" must be integer greater than 0');
                    return;
                end

                % validate n_data is number
                
                data.n = n_data;
                
                % validate matrixes # sum( sum( isnan(M) ) ) == 0
                data.method = 'mean_based';
                if strcmp(selected_model, 'User Defined')
                    data.method = 'user defined';
                end
                data.parameters = user_defined_params;
                    
                delete(gen_dialog);
            end
            
            function cancel_action(src, ev)
                delete(gen_dialog);
            end
            
            mean_based.Callback = @mean_based_action;
            user_defined.Callback = @user_defined_action;
            help_button.Callback = @help_action;
            generate_button.Callback = @generate_action;
            cancel_button.Callback = @cancel_action;

            uiwait(gen_dialog);
        end
        
        function data = get_select_vars(obj, btn_txt)
            data = {};
            checks = {};

            function select(src, ~)
                selected = fields(checks);
                for i=1:size(selected,1)
                    check = getfield(checks, selected{i});
                    check.Value = 1;
                end
            end
            
            select_dialog = figure(...
                'units', 'pixels',...
                'Position', [100 100 250 10],...
                'Toolbar', 'none',...
                'Menubar', 'none',...
                'Name', 'Variables',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            y = 40;
            % opt values
            if ~isempty(obj.study.opt)
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Training options (opt)',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['opt']);
                checks = setfield(checks, 'opt', aux);
                y = y + 30;
            end
            % viterbi sequence
            if ~isempty(obj.study.viterbi_results)
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Viterbi Sequence',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['viterbi']);
                checks = setfield(checks, 'viterbi', aux);
                y = y + 30;
            end
            % generated data
            if ~isempty(obj.study.generated_data)
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Generated Data',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['generated']);
                checks = setfield(checks, 'generated', aux);
                y = y + 30;
            end
            % decoded data
            if ~isempty(obj.study.decode_results)
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Decode Results',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['decode']);
                checks = setfield(checks, 'decode', aux);
                y = y + 30;
            end
            % results
            if ~isempty(obj.study.training_results)
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Training Results',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['results']);
                checks = setfield(checks, 'results', aux);
                y = y + 30;
            end
            % model
            if ~isempty(obj.study.model)
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Model',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['model']);
                checks = setfield(checks, 'model', aux);
                y = y + 30;
            end
            % training data
            if ~isempty(obj.study.data)
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Loaded Data',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['data']);
                checks = setfield(checks, 'data', aux);
                y = y + 30;
            end

            uicontrol(select_dialog,...
                'Style', 'pushbutton',...
                'String', 'Select All',...
                'Position', [45 y 160 20],...
                'Callback', @select);
            
            select_dialog.Position(4) = y+30;
            
            % ok and cancel buttons
            ok_button = uicontrol(select_dialog,...
                'Style', 'pushbutton',...
                'String', btn_txt,...
                'ForegroundColor', [0.0 0.8 0.0],...
                'Position', [130 10 50 20]);
            cancel_button = uicontrol(select_dialog,...
                'Style', 'pushbutton',...
                'String', 'Cancel',...
                'ForegroundColor', [0.8 0.0 0.0],...
                'Position', [190 10 50 20]);

            function ok_action(~, ~)
                aux = fields(checks);
                for i=1:size(aux,1)
                    check = getfield(checks, aux{i});
                    data = setfield(data, aux{i}, check.Value);
                end
                delete(select_dialog);
            end
            
            function cancel_action(~, ~)
                delete(select_dialog);
            end
            
            ok_button.Callback = @ok_action;
            cancel_button.Callback = @cancel_action;

            uiwait(select_dialog);
        end
        
        function data = get_select_vars2(~, btn_txt, list)
            data = {};
            checks = {};

            function select(src, ~)
                selected = fields(checks);
                for i=1:size(selected,1)
                    check = getfield(checks, selected{i});
                    check.Value = 1;
                end
            end
            
            select_dialog = figure(...
                'units', 'pixels',...
                'Position', [100 100 250 10],...
                'Toolbar', 'none',...
                'Menubar', 'none',...
                'Name', 'Variables',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            y = 40;
            % opt values
            if any(strcmp(list, 'opt'))
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Trainig options (opt)',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['opt']);
                checks = setfield(checks, 'opt', aux);
                y = y + 30;
            end
            % viterbi sequence
            if any(strcmp(list, 'viterbi_results'))
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Viterbi Sequence',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['viterbi']);
                checks = setfield(checks, 'viterbi', aux);
                y = y + 30;
            end
            % generated data
            if any(strcmp(list, 'generated_data'))
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Generated Data',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['generated']);
                checks = setfield(checks, 'generated', aux);
                y = y + 30;
            end
            % decoded data
            if any(strcmp(list, 'decode_results'))
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Decode Results',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['decode']);
                checks = setfield(checks, 'decode', aux);
                y = y + 30;
            end
            % results
            if any(strcmp(list, 'training_results'))
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Training Results',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['results']);
                checks = setfield(checks, 'results', aux);
                y = y + 30;
            end
            % model
            if any(strcmp(list, 'model'))
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Model',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['model']);
                checks = setfield(checks, 'model', aux);
                y = y + 30;
            end
            % training data
            if any(strcmp(list, 'data'))
                aux = uicontrol(select_dialog,...
                    'Style', 'checkbox',...
                    'String', 'Loaded Data',...
                    'Position', [25 y 210 20],...
                    'Value', 0,...
                    'UserData', ['data']);
                checks = setfield(checks, 'data', aux);
                y = y + 30;
            end

            uicontrol(select_dialog,...
                'Style', 'pushbutton',...
                'String', 'Select All',...
                'Position', [45 y 160 20],...
                'Callback', @select);
            
            select_dialog.Position(4) = y+30;
            
            % ok and cancel buttons
            ok_button = uicontrol(select_dialog,...
                'Style', 'pushbutton',...
                'String', btn_txt,...
                'ForegroundColor', [0.0 0.8 0.0],...
                'Position', [130 10 50 20]);
            cancel_button = uicontrol(select_dialog,...
                'Style', 'pushbutton',...
                'String', 'Cancel',...
                'ForegroundColor', [0.8 0.0 0.0],...
                'Position', [190 10 50 20]);

            function ok_action(~, ~)
                aux = fields(checks);
                for i=1:size(aux,1)
                    check = getfield(checks, aux{i});
                    data = setfield(data, aux{i}, check.Value);
                end
                delete(select_dialog);
            end
            
            function cancel_action(~, ~)
                delete(select_dialog);
            end
            
            ok_button.Callback = @ok_action;
            cancel_button.Callback = @cancel_action;

            uiwait(select_dialog);
        end
        
        % check exists opt and ask to create it
        function data = check_opt(obj)
            data = 0;
            
            if isempty(obj.study.opt)
                choice = questdlg({'"opt" variable does not exist,', 'Do you want to create it by default?'},...
                    'Missing options',...
                    'Yes', 'No', 'Yes');
                if strcmp(choice, 'Yes')
                    obj.load_opt();
                    data = 1;
                end
            else
                data = 1;
                return;
            end
        end
        
        % emission and duration distribution panel creation
        function g = get_emission_duration_panel(~, parent, name, model, is_prior)
            function modify_cell(src, ev, is_prior, cbox)
                is_editable = is_prior && (cbox.Value>1);

                % nobody generate this event so skip
                if isempty(ev.Indices)
                    return
                end
                
                col = ev.Indices(2); % n_states
                row = ev.Indices(1); % parameter name
                
                text = 'weiv';
                if is_editable
                    text = 'tide';
                end
                [aux{1:size(src.Data,1),1:size(src.Data,2)}] = deal(text);
                src.Data = aux;
                text = 'view';
                if is_editable
                    text = 'edit';
                end
                [aux{1:size(src.Data,1),1:size(src.Data,2)}] = deal(text);
                src.Data = aux;

                row_name = src.RowName{row};
                the_model = src.UserData;
                sub_dialog = figure(...
                    'units', 'pixels',...
                    'Position', [100 100 182 166],...
                    'Toolbar', 'none',...
                    'Menubar', 'none',...
                    'Name', '',...
                    'numbertitle', 'off',...
                    'resize', 'off');
               
                f1 = getfield(the_model{col}, row_name);
                cw = num2cell(40*ones(1, size(f1, 2)));
                ed = false(1, size(f1, 2));
                if is_editable
                    ed = true(1, size(f1, 2));
                end

                my_table = uitable(sub_dialog,...
                    'ColumnEditable', ed,...
                    'Data', f1,...
                    'ColumnWidth', cw,...
                    'Position', [10 40 162 86]);

                text = '';
                if ~is_prior
                    text = ['State N°' num2str(col) ' / '];
                end
                text = [text 'Parameter: ' row_name];
                the_label = uicontrol(sub_dialog,...
                    'Position', [10 136 162 20],...
                    'Style', 'text',...
                    'HorizontalAlignment', 'left',...
                    'String', text);

                function load_action(~, ~)
                    % retrieve var from global workspace
                    vars = evalin('base', 'who()');
                    % ask the user to select one
                    [indx,tf] = listdlg('PromptString', 'Select Data Variable:',...
                        'SelectionMode', 'single',...
                        'ListString', vars);
                    if ~tf
                        % cancel the loading by close or press cancel
                        return;
                    end
                    % get the value
                    value = evalin('base', vars{indx});
                    my_table.Data = value;
                end
                    
                function ok_action(~, ~)
                    isn = find(isnan(my_table.Data));
                    if ~isempty(isn) && size(isn, 1) > 0
                        obj.my_warning('Matrix (or vector) cannot contains NaN values');
                        return;
                    end
                    the_model{col} = setfield(the_model{col}, row_name, my_table.Data);
                    src.UserData = the_model;
                    delete(sub_dialog);
                end
                
                function cancel_action(~, ~)
                    delete(sub_dialog);
                end

                if is_editable
                    load_button = uicontrol(sub_dialog,...
                        'Position', [122 136 50 20],...
                        'Style', 'pushbutton',...
                        'String', 'load');

                    load_button.Callback = @load_action;
                    
                    the_label.Position(2) = the_label.Position(2)+30;
                    sub_dialog.Position(4) = sub_dialog.Position(4)+30;
                    
                    ok_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'OK',...
                        'ForegroundColor', [0.0 0.8 0.0],...
                        'Position', [62 10 50 20]);
                    cancel_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'Cancel',...
                        'ForegroundColor', [0.8 0.0 0.0],...
                        'Position', [122 10 50 20]);
                    ok_button.Callback = @ok_action;
                    cancel_button.Callback = @cancel_action;
                else
                    ok_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'OK',...
                        'ForegroundColor', [0.0 0.8 0.0],...
                        'Position', [122 10 50 20]);
                    ok_button.Callback = @cancel_action;
                end
            end
            
            % model name
            m_name = class(model);
            n = size(m_name, 2);
            i = strfind(m_name, '.');
            m_name = m_name(i+1:n);
            % distribution
            g.ds_panel = uipanel(parent,...
                'units', 'pixels',...
                'Title', [name ' Model'],...
                'FontWeight', 'bold',...
                'Position', [0 0 500 135]);
            distribution_model = model.prior;
            if ~is_prior
                distribution_model = model.posterior;
            end
            
            g.ds_name = uicontrol(g.ds_panel,...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'Position', [10 95 250 20],...
                'String', m_name);
            x = 10;
            if is_prior
                g.ds_combo = uicontrol(g.ds_panel,...
                    'Style', 'popup',...
                    'HorizontalAlignment', 'left',...
                    'Position', [x 70 100 25],...
                    'String', {'Non Informative', 'User Defined'});
                x = x + 110;
            else
                g.ds_combo = {};
            end

            % fields of the model
            params = fields(distribution_model);
            % visit every field
            g.ds_params = {};
            for i=1:size(params,1)
                key = params{i};
                value = getfield(distribution_model, key);
                my_panel = uipanel(g.ds_panel,...
                    'Title', key,...
                    'units', 'pixels',...
                    'Position', [x 10 180 90]);
                
                % number of states in the model
                n_states = size(value, 2);
                % parameters of states (same for everyone)
                state_params = fields(value{1});
                % matrix of edit string
                [model_data{1:size(state_params,1),1:n_states}] = deal('view');
                % width for columns
                [cw{1:n_states}] = deal(30);
                % finally the table
                my_table = uitable(my_panel,...
                    'RowName', state_params,...
                    'Data', model_data,...
                    'ColumnWidth', cw,...
                    'Position', [10 10 160 60],...
                    'CellSelectionCallback', {@modify_cell, is_prior, g.ds_combo},...
                    'UserData', value);

                g.ds_params = setfield(g.ds_params, key, my_table);
                % move pointer
                x = x+190;
                clear('model_data');
            end
            
            g.ds_panel.Position(3) = x;
        end
        
        % emission and duration distribution panel creation
        function g = get_transition_initial_panel(~, parent, name, model, is_prior)
            function modify_cell(src, ev, is_prior, cbox)
                is_editable = is_prior && (cbox.Value>1);

                % nobody generate this event so skip
                if isempty(ev.Indices)
                    return
                end
                
                col = ev.Indices(2); % n_states
                row = ev.Indices(1); % parameter name
                
                text = 'weiv';
                if is_editable
                    text = 'tide';
                end
                [aux{1:size(src.Data,1),1:size(src.Data,2)}] = deal(text);
                src.Data = aux;
                text = 'view';
                if is_editable
                    text = 'edit';
                end
                [aux{1:size(src.Data,1),1:size(src.Data,2)}] = deal(text);
                src.Data = aux;

                row_name = src.RowName{row};
                the_model = src.UserData;
                sub_dialog = figure(...
                    'units', 'pixels',...
                    'Position', [100 100 182 166],...
                    'Toolbar', 'none',...
                    'Menubar', 'none',...
                    'Name', '',...
                    'numbertitle', 'off',...
                    'resize', 'off');
               
                f1 = getfield(the_model, row_name);
                cw = num2cell(40*ones(1, size(f1, 2)));
                ed = false(1, size(f1, 2));
                if is_editable
                    ed = true(1, size(f1, 2));
                end

                my_table = uitable(sub_dialog,...
                    'ColumnEditable', ed,...
                    'Data', f1,...
                    'ColumnWidth', cw,...
                    'Position', [10 40 162 86]);

                the_label = uicontrol(sub_dialog,...
                    'Position', [10 136 162 20],...
                    'Style', 'text',...
                    'HorizontalAlignment', 'left',...
                    'String', ['Parameter: ' row_name]);

                function load_action(~, ~)
                    % retrieve var from global workspace
                    vars = evalin('base', 'who()');
                    % ask the user to select one
                    [indx,tf] = listdlg('PromptString', 'Select Data Variable:',...
                        'SelectionMode', 'single',...
                        'ListString', vars);
                    if ~tf
                        % cancel the loading by close or press cancel
                        return;
                    end
                    % get the value
                    value = evalin('base', vars{indx});
                    my_table.Data = value;
                end
                    
                function ok_action(~, ~)
                    isn = find(isnan(my_table.Data));
                    if ~isempty(isn) && size(isn, 1) > 0
                        obj.my_warning('Matrix (or vector) cannot contains NaN values');
                        return;
                    end
                    the_model = setfield(the_model, row_name, my_table.Data);
                    src.UserData = the_model;
                    delete(sub_dialog);
                end
                
                function cancel_action(~, ~)
                    delete(sub_dialog);
                end

                if is_editable
                    load_button = uicontrol(sub_dialog,...
                        'Position', [122 136 50 20],...
                        'Style', 'pushbutton',...
                        'String', 'load');

                    load_button.Callback = @load_action;
                    
                    the_label.Position(2) = the_label.Position(2)+30;
                    sub_dialog.Position(4) = sub_dialog.Position(4)+30;
                    
                    ok_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'OK',...
                        'ForegroundColor', [0.0 0.8 0.0],...
                        'Position', [62 10 50 20]);
                    cancel_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'Cancel',...
                        'ForegroundColor', [0.8 0.0 0.0],...
                        'Position', [122 10 50 20]);
                    ok_button.Callback = @ok_action;
                    cancel_button.Callback = @cancel_action;
                else
                    ok_button = uicontrol(sub_dialog,...
                        'Style', 'pushbutton',...
                        'String', 'OK',...
                        'ForegroundColor', [0.0 0.8 0.0],...
                        'Position', [122 10 50 20]);
                    ok_button.Callback = @cancel_action;
                end
            end
            
            % model name
            m_name = class(model);
            n = size(m_name, 2);
            i = strfind(m_name, '.');
            m_name = m_name(i+1:n);
            % distribution
            g.ds_panel = uipanel(parent,...
                'units', 'pixels',...
                'Title', [name ' Model'],...
                'FontWeight', 'bold',...
                'Position', [0 0 500 135]);
            if ~model.ndim
                g.ds_panel.ForegroundColor = [0.7 0.7 0.7];
            else
                g.ds_panel.ForegroundColor = [0 0 0];
            end

            distribution_model = model.prior;
            if ~is_prior
                distribution_model = model.posterior;
            end
            g.ds_name = uicontrol(g.ds_panel,...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'Position', [10 95 250 20],...
                'String', m_name);
            if ~model.ndim
                g.ds_name.Enable = 'off';
            end
            x = 10;
            if is_prior
                g.ds_combo = uicontrol(g.ds_panel,...
                    'Style', 'popup',...
                    'HorizontalAlignment', 'left',...
                    'Position', [x 70 100 25],...
                    'String', {'Non Informative', 'User Defined'});
                if ~model.ndim
                    g.ds_combo.Enable = 'off';
                end
                x = x + 110;
            else
                g.ds_combo = {};
            end

            % fields of the model
            params = fields(distribution_model);
            % visit every field
            g.ds_params = {};
            for i=1:size(params,1)
                key = params{i};
                value = getfield(distribution_model, key);
                my_panel = uipanel(g.ds_panel,...
                    'Title', key,...
                    'units', 'pixels',...
                    'Position', [x 10 180 90]);
                if ~model.ndim
                    my_panel.ForegroundColor = [0.7 0.7 0.7];
                else
                    my_panel.ForegroundColor = [0 0 0];
                end
                
                % parameters of states (same for everyone)
                state_params = fields(value);
                % matrix of edit string
                [model_data{1:size(state_params,1),1}] = deal('view');%
                % width for columns
                cw = {30};
                % finally the table
                my_table = uitable(my_panel,...
                    'RowName', state_params,...
                    'Data', model_data,...
                    'ColumnWidth', cw,...
                    'Position', [10 10 160 60],...
                    'CellSelectionCallback', {@modify_cell, is_prior, g.ds_combo},...
                    'UserData', value);
                if ~model.ndim
                    my_table.Enable = 'off';
                else
                    my_table.CellSelectionCallback = {@modify_cell, is_prior, g.ds_combo};
                end
                g.ds_params = setfield(g.ds_params, key, my_table);
                % move pointer
                x = x+190;
                clear('model_data');
            end
            
            g.ds_panel.Position(3) = x;
        end
        
        % retrieve the route for the eeg hardware coordinates
        function g = select_eeg_coords(obj)
            g = {};
            path = '';
            
            % dialog
            eeg_dialog = figure(...
                'units', 'pixels',...
                'Position', [100 100 300 100],...
                'Toolbar', 'none',...
                'Menubar', 'none',...
                'Name', 'Select EEG Coordinates',...
                'numbertitle', 'off',...
                'resize', 'off');
            
            function get_eeg_file(src, ev)
                if src == browse_radio
                    browse_btn.Enable = 'on';
                else
                    browse_btn.Enable = 'off';
                end
                path = src.String;
                if strcmp(path, '<Seleccionar archivo>')
                    path = '';
                end
            end
            
            function browse_action(src, ev)
                [selected_file, path] = uigetfile({'*.xyz',  'XYZ Files (*.xyz)'});
                if ~isempty(selected_file)
                    path = [path selected_file];
                    browse_radio.String = ['[...]\' selected_file];
                end
                disp('browsing');
            end
            
            % button group for radios
            eeg_options = uibuttongroup(eeg_dialog,...
                'Title', 'Coordinates File',...
                'units', 'pixels',...
                'Position', [10 40 280 80]);
            browse_radio = uicontrol(eeg_options,...
                'Style', 'radiobutton',...
                'Position', [10 10 200 20],...
                'String', '<Seleccionar archivo>',...
                'Callback', @get_eeg_file);
            % browse radio
            browse_btn = uicontrol(eeg_options,...
                'Style', 'pushbutton',...
                'Position', [210 10 60 20],...
                'String', 'Browse',...
                'Callback', @browse_action);
            % default elements
            flipped = flip(obj.default_eeg);
            for i=1:size(obj.default_eeg,2)
                uicontrol(eeg_options,...
                    'Style', 'radiobutton',...
                    'Position', [10 (10+25*i) 260 20],...
                    'String', flipped{i},...
                    'Callback', @get_eeg_file);
            end
            eeg_options.Position(4) = 45+25*size(obj.default_eeg,2);
            eeg_dialog.Position(4) = eeg_options.Position(4)+40;
            
            function ok_action(~, ~)
                g = path;
                delete(eeg_dialog);
            end
            
            function cancel_action(~, ~)
                g = {};
                delete(eeg_dialog);
            end
            
            ok_button = uicontrol(eeg_dialog,...
                'Style', 'pushbutton',...
                'String', 'OK',...
                'ForegroundColor', [0.0 0.8 0.0],...
                'Position', [180 10 50 20],...
                'Callback', @ok_action);
            cancel_button = uicontrol(eeg_dialog,...
                'Style', 'pushbutton',...
                'String', 'Cancel',...
                'ForegroundColor', [0.8 0.0 0.0],...
                'Position', [240 10 50 20],...
                'Callback', @cancel_action);
            
            uiwait(eeg_dialog);
        end
        
        % refresh display of study variables
        function show_study(obj)
            % if there is no data at all
            at_least_one = 0;
            % DATA
            if ~isempty(obj.study.data)
                at_least_one = 1;
                % panel color and text
                obj.gui.data_panel.ForegroundColor = [0.0 0.8 0.0];
                obj.gui.data_panel.Title = ['Name : [ ' obj.study.metadata.name ' ]'];
                % labels
                obj.gui.data_filename_label.ForegroundColor = [0.0 0.0 0.0];
                obj.gui.data_filename_value.ForegroundColor = [0.0 0.0 0.8];
                pieces = strsplit(obj.study.metadata.filename, '\');
                obj.gui.data_filename_value.String = pieces(size(pieces,2));
                obj.gui.data_dimension_label.ForegroundColor = [0.0 0.0 0.0];
                obj.gui.data_dimension_value.ForegroundColor = [0.0 0.0 0.8];
                obj.gui.data_dimension_value.String = num2str(obj.study.metadata.dimension);
                obj.gui.data_length_value.ForegroundColor = [0.0 0.0 0.8];
                obj.gui.data_length_label.ForegroundColor = [0.0 0.0 0.0];
                if obj.study.metadata.multisubject;
                    obj.gui.data_length_label.Visible = 'off';
                    obj.gui.data_length_value.String = 'MultiSubject';
                else
                    obj.gui.data_length_value.String = num2str(obj.study.metadata.length);
                end
            else
                % panel color and text
                obj.gui.data_panel.ForegroundColor = [0.8 0.0 0.0];
                obj.gui.data_panel.Title = 'Data not Loaded';
                % labels
                obj.gui.data_filename_label.ForegroundColor = [0.5 0.5 0.5];
                obj.gui.data_filename_value.String = '';
                obj.gui.data_dimension_label.ForegroundColor = [0.5 0.5 0.5];
                obj.gui.data_dimension_value.String = '';
                obj.gui.data_length_label.ForegroundColor = [0.5 0.5 0.5];
                obj.gui.data_length_label.Visible = 'on';
                obj.gui.data_length_value.String = '';
            end
            % MODEL
            if ~isempty(obj.study.model)
                if strcmp(class(obj.study.model), 'multimodel')
                    my_model = obj.study.model.matrixmodel(1,1,1);
                else
                    my_model = obj.study.model;
                end
                at_least_one = 1;
                % panel color and text
                obj.gui.model_panel.ForegroundColor = [0.0 0.8 0.0];
                obj.gui.model_panel.Title = ['Model [ ' upper(class(obj.study.model)) ' ]'];
                % labels
                obj.gui.prior_label.ForegroundColor = [0.0 0.0 0.0];
                obj.gui.prior_value.ForegroundColor = [0.0 0.8 0.0];
                obj.gui.prior_value.String = 'CREATED';
                obj.gui.posterior_label.ForegroundColor = [0.0 0.0 0.0];
                if my_model.emis_model.posteriorfull
                    obj.gui.posterior_value.ForegroundColor = [0.0 0.8 0.0];
                    obj.gui.posterior_value.String = 'CREATED';
                else
                    obj.gui.posterior_value.ForegroundColor = [0.8 0.0 0.0];
                    obj.gui.posterior_value.String = 'NOT CREATED';
                end
                % # states
                obj.gui.nstates_label.ForegroundColor = [0.0 0.0 0.0];
                if my_model.nstates > 0
                    obj.gui.nstates_value.ForegroundColor = [0.0 0.8 0.0];
                    obj.gui.nstates_value.String = num2str(my_model.nstates);
                else
                    obj.gui.nstates_value.ForegroundColor = [0.8 0.0 0.0];
                    obj.gui.nstates_value.String = 'undefined';
                end
            else
                % panel color and text
                obj.gui.model_panel.ForegroundColor = [0.8 0.0 0.0];
                obj.gui.model_panel.Title = 'Model not Created';
                % labels
                obj.gui.prior_label.ForegroundColor = [0.5 0.5 0.5];
                obj.gui.prior_value.String = '';
                obj.gui.posterior_label.ForegroundColor = [0.5 0.5 0.5];
                obj.gui.posterior_value.String = '';
                obj.gui.nstates_label.ForegroundColor = [0.5 0.5 0.5];
                obj.gui.nstates_value.String = '';
            end
            
            if ~at_least_one
                obj.gui.tab_status.ForegroundColor = [0.8 0.0 0.0];
            else
                obj.gui.tab_status.ForegroundColor = [0.0 0.8 0.0];
            end
        end
        
        function classes_list=classes_ls(~, path)
            complete_path = ['+' path];
            elements = ls(complete_path);
            classes_list = {};
            for e=3:size(elements, 1)
                a = strtrim(elements(e,:));
                if strcmp(a(size(a,2)-1:size(a,2)), '.m')
                    b = a(1:size(a,2)-2);
                    classes_list = [classes_list b];
                end
            end
        end
    
        function set_status_bar(obj, msj)
            obj.gui.status_bar.String = msj;
            drawnow();
        end

        function menu_enable_state(obj)
            ti_multimod = 0;
            if strcmp(class(obj.study.model), 'multimodel')
                my_model = obj.study.model.matrixmodel(1,1,1);
                ti_multimod = 1;
            else
                my_model = obj.study.model;
            end
            ti_data = ~isempty(obj.study.data);
            ti_mod = ~isempty(my_model);
            ti_posterior = 0;
            if ti_mod
                ti_posterior = my_model.emis_model.posteriorfull();
            end
            ti_has_dur = strcmp(class(my_model), 'hsmm');
            ti_train = ~isempty(obj.study.training_results);
            ti_gen = ~isempty(obj.study.generated_data);
            ti_seq = ~isempty(obj.study.viterbi_results);
            ti_dec = ~isempty(obj.study.decode_results);
            ti_ws = (ti_data || ti_mod || ti_train || ti_gen || ti_seq || ti_dec);
            
            % export data
            if ti_gen
                obj.gui.export_menu.Enable = 'on';
            else
                obj.gui.export_menu.Enable = 'off';
            end
            % save and clear study
            if ti_ws
                obj.gui.save_study_menu.Enable = 'on';
                obj.gui.save_as_study_menu.Enable = 'on';
                obj.gui.clear_study_menu.Enable = 'on';
            else
                obj.gui.save_study_menu.Enable = 'off';
                obj.gui.save_as_study_menu.Enable = 'off';
                obj.gui.clear_study_menu.Enable = 'off';
            end            
            % prior
            if ti_mod
                obj.gui.prior_dist_menu.Enable = 'on';
                obj.gui.plot_gamma_menu.Enable = 'on';
            else
                obj.gui.prior_dist_menu.Enable = 'off';
                obj.gui.plot_gamma_menu.Enable = 'off';
            end            
            % posterior, estimate parameters
            if ti_posterior
                obj.gui.posterior_dist_menu.Enable = 'on';
            else
                obj.gui.posterior_dist_menu.Enable = 'off';
            end
            % train model
            if ti_mod && ti_data && ~ti_multimod
                obj.gui.train_model_menu.Enable = 'on';
            else
                obj.gui.train_model_menu.Enable = 'off';
            end
            % generate
            if ti_mod
                obj.gui.generate_data_menu.Enable = 'on';
            else
                obj.gui.generate_data_menu.Enable = 'off';
            end            
            % viterbi, estimate parameters
            if ti_posterior && ti_data
                obj.gui.sequence_estimation_menu.Enable = 'on';
                obj.gui.estimate_parameters_menu.Enable = 'on';
            else
                obj.gui.sequence_estimation_menu.Enable = 'off';
                obj.gui.estimate_parameters_menu.Enable = 'off';
            end            

            % decode
            if ti_train && ti_data
                obj.gui.data_decode_menu.Enable = 'on';
            else
                obj.gui.data_decode_menu.Enable = 'of';
            end            
            % sequence and ocupancy plot
            if ti_seq
                obj.gui.plot_sequence_menu.Enable = 'on';
                obj.gui.plot_occupancy_menu.Enable = 'on';
            else
                obj.gui.plot_sequence_menu.Enable = 'off';
                obj.gui.plot_occupancy_menu.Enable = 'off';
            end
            % duration histogram plot
            if ti_seq && ti_has_dur
                obj.gui.plot_histogram_menu.Enable = 'on';
            else
                obj.gui.plot_histogram_menu.Enable = 'off';
            end
            % topography plot
            if ti_posterior
                obj.gui.plot_topo_menu.Enable = 'on';
            else
                obj.gui.plot_topo_menu.Enable = 'off';
            end
            % emission covariance plot
            if ti_posterior
                obj.gui.plot_posterior_dur_menu.Enable = 'on';
                obj.gui.plot_cov_matrix_menu.Enable = 'on';
            else
                obj.gui.plot_posterior_dur_menu.Enable = 'off';
                obj.gui.plot_cov_matrix_menu.Enable = 'off';
            end
            % free energy
            if ti_train
                obj.gui.plot_free_energy_menu.Enable = 'on';
            else
                obj.gui.plot_free_energy_menu.Enable = 'off';
            end
            
            % emission covariance plot
            if ti_gen
                obj.gui.plot_gen_sequence_menu.Enable = 'on';
                obj.gui.plot_gen_occupancy_menu.Enable = 'on';
                obj.gui.plot_gen_durhist_menu.Enable = 'on';
                obj.gui.plot_gen_topo_menu.Enable = 'on';
                obj.gui.plot_gen_covmat_menu.Enable = 'on';
            else
                obj.gui.plot_gen_sequence_menu.Enable = 'off';
                obj.gui.plot_gen_occupancy_menu.Enable = 'off';
                obj.gui.plot_gen_durhist_menu.Enable = 'off';
                obj.gui.plot_gen_topo_menu.Enable = 'off';
                obj.gui.plot_gen_covmat_menu.Enable = 'off';
            end
            
        end
        
        function multisubject_selector(obj, panel, cb)
            my_panel = uipanel(panel,...
                'units', 'pixels',...
                'Title', ['Model Selector'],...
                'FontWeight', 'bold',...
                'Position', [10 panel.Position(4) panel.Position(3)-20 50]);
            panel.Position(4) = panel.Position(4)+60;
            
            uicontrol(my_panel,...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'Position', [10 5 30 20],...
                'String', 'cond');
            combo_cond = uicontrol(my_panel,...
                'Style', 'popup',...
                'HorizontalAlignment', 'left',...
                'Position', [50 10 40 20],...
                'String', 1:size(obj.study.data.cond, 2));
            uicontrol(my_panel,...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'Position', [100 5 30 20],...
                'String', 'subj');
            combo_subj = uicontrol(my_panel,...
                'Style', 'popup',...
                'HorizontalAlignment', 'left',...
                'Position', [140 10 40 20],...
                'String', 1:size(obj.study.data.cond(1).subj, 2));
            uicontrol(my_panel,...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'Position', [190 5 30 20],...
                'String', 'block');
            combo_block = uicontrol(my_panel,...
                'Style', 'popup',...
                'HorizontalAlignment', 'left',...
                'Position', [230 10 40 20],...
                'String', 1:size(obj.study.data.cond(1).subj(1).block, 2));

            function my_callback()
                v_cond = combo_cond.Value;
                v_subj = combo_subj.Value;
                v_block = combo_block.Value;
                cb(v_cond, v_subj, v_block);
            end
                
            function cb_cond(~, ~)
                v_cond = combo_cond.Value;
                combo_subj.Value = 1;
                combo_subj.String = 1:size(obj.study.data.cond(v_cond).subj, 2);
                combo_block.Value = 1;
                combo_block.String = 1:size(obj.study.data.cond(v_cond).subj(1).block, 2);
                my_callback();
            end
            
            function cb_subj(~, ~)
                v_cond = combo_cond.Value;
                v_subj = combo_subj.Value;
                combo_block.Value = 1;
                combo_block.String = 1:size(obj.study.data.cond(v_cond).subj(v_subj).block, 2);
                my_callback();
            end
            
            function cb_block(~, ~)
                my_callback();
            end
            
            combo_cond.Callback = @cb_cond;
            combo_subj.Callback = @cb_subj;
            combo_block.Callback = @cb_block;
        end
        
        function my_error(obj, err)
            errordlg({'Ooops! Something happened...';'Watch the console for further information'},'Error');
            disp('########################');
            disp(['# ' datestr(datetime('now')) ' #']);
            disp('########################');disp(' ');
            disp(err.getReport());
            obj.set_status_bar('Something went wrong! =(')
        end
        
        function my_warning(~, msj)
            warndlg(msj,'Warning');
        end
    
        function update_title(obj)
            %Cambio David "BDS" a "BSD"
            txt = 'BSD [';
            disp(obj.ws_filename);
            if isempty(obj.ws_filename)
                txt = [txt 'untitled'];
            else
                txt = [txt obj.ws_filename];
            end
            if ~obj.is_saved
                txt = [txt '*'];
            end
            txt = [txt ']'];
            
            obj.gui.fig.Name = txt;
        end
        
        function set_saved_status(obj, state)
            obj.show_study();
            obj.is_saved = state;
            obj.menu_enable_state();
            obj.update_title();
        end
        
        function result=is_save_query(obj, title, question)
            result=1;
            % if there is something
            if ~obj.is_saved
                choice = questdlg(question,...
                title,...
                'Yes', 'No', 'Yes');
                if isempty(choice)
                    result=0;
                    return;
                end
                if strcmp(choice, 'Yes')
                    result = obj.save_study_action(0, 0, 0);
                end
            end
        end
        
        function leave(obj, ~, ~)
            if obj.is_save_query('Not saved...', 'Do you want to save the changes?')
                delete(obj.gui.fig);
            end
        end
    end
    
end