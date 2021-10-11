classdef retroData
    
    % Data and parameter class for retrospective app
    
    properties
        
        % Data
        data
        mrdFooter
        rprFile
        filename
        includeWindow                                                  
        excludeWindow 
        nrKlines
        
        % Sequence parameters
        PPL
        NO_SAMPLES = 1
        NO_VIEWS = 1
        NO_VIEWS_2 = 1
        EXPERIMENT_ARRAY = 1
        nr_repetitions = 1
        NO_AVERAGES = 1
        NO_SLICES = 1
        SLICE_THICKNESS = 1
        SLICE_SEPARATION = 1
        SLICE_INTERLEAVE = 1
        r_angle_var
        p_angle_var
        s_angle_var
        nr_coils = 1
        FOV = 30
        PHASE_ORIENTATION = 0
        FOVf = 8
        aspectratio = 1
        alpha = 20
        te = 2
        te_us = 0
        TE
        tr = 10
        tr_extra_us = 0
        TR
        ti = 1000
        VFA_angles
        VFA_size = 0
        frame_loop_on
        radial_on = 0
        slice_nav = 0
        date
        pixelshift1 = 0
        pixelshift2 = 0
        coil_scaling = 1
        scanner = 'MRS'
        
        % K-space trajectory related
        pe1_order = 3
        pe2_centric_on = 1
        pe2_traj = 0
        gp_var_mul
        
        % Navigator related
        no_samples_nav = 10
        primary_navigator_point = 10
        nr_nav_points_discarded = 35
        nr_nav_points_used = 5
        
        % Final dimensions
        dimx
        dimy
        dimz
        nrKsteps
        
        % Flags
        rpr_flag = false;
        validData_flag = false;
        multiCoil_flag = false
        multi2D_flag = false
        vfaData_flag = false
        
        % Data and reconstruction type
        dataType = '2D'
        reco_guess = 'systolic function'
        
    end
    
    
    
    methods
        
        % ---------------------------------------------------------------------------------
        % Object constructor
        % ---------------------------------------------------------------------------------
        function obj = retroData(parameter,mridata)
            
            if nargin == 2
                
                obj.data = mridata;
                
                if isfield(parameter,'filename')
                    obj.filename = parameter.filename;
                end
                
                if isfield(parameter,'PPL')
                    obj.PPL = parameter.PPL;
                end
                
                if isfield(parameter,'NO_SAMPLES')
                    obj.NO_SAMPLES = parameter.NO_SAMPLES;
                end
                
                if isfield(parameter,'NO_VIEWS')
                    obj.NO_VIEWS = parameter.NO_VIEWS;
                end
                
                if isfield(parameter,'NO_VIEWS_2')
                    obj.NO_VIEWS_2 = parameter.NO_VIEWS_2;
                end
                
                if isfield(parameter,'EXPERIMENT_ARRAY')
                    obj.EXPERIMENT_ARRAY = parameter.EXPERIMENT_ARRAY;
                    obj.nr_repetitions = parameter.EXPERIMENT_ARRAY;
                end
                
                if isfield(parameter,'NO_AVERAGES')
                    obj.NO_AVERAGES = parameter.NO_AVERAGES;
                end
                
                if isfield(parameter,'NO_SLICES')
                    obj.NO_SLICES = parameter.NO_SLICES;
                end
                
                if isfield(parameter,'SLICE_THICKNESS')
                    obj.SLICE_THICKNESS= parameter.SLICE_THICKNESS;
                end
                
                if isfield(parameter,'SLICE_SEPARATION')
                    obj.SLICE_SEPARATION = parameter.SLICE_SEPARATION;
                end
                
                if isfield(parameter,'SLICE_INTERLEAVE')
                    obj.SLICE_INTERLEAVE = parameter.SLICE_INTERLEAVE;
                end
                
                if isfield(parameter,'r_angle_var')
                    obj.r_angle_var = parameter.r_angle_var;
                end
                
                if isfield(parameter,'p_angle_var')
                    obj.p_angle_var = parameter.p_angle_var;
                end
                
                if isfield(parameter,'s_angle_var')
                    obj.s_angle_var = parameter.s_angle_var;
                end
                
                if isfield(parameter,'nr_coils')
                    obj.nr_coils = parameter.nr_coils;
                    if obj.nr_coils > 1
                        obj.multiCoil_flag = true;
                    else
                        obj.multiCoil_flag = false;
                    end
                end
                
                if isfield(parameter,'FOV')
                    obj.FOV = parameter.FOV;
                end
                
                if isfield(parameter,'PHASE_ORIENTATION')
                    obj.PHASE_ORIENTATION = parameter.PHASE_ORIENTATION;
                end
                
                if isfield(parameter,'FOVf')
                    obj.FOVf = parameter.FOVf;
                end
                
                obj.aspectratio = obj.FOVf/8;
                
                if isfield(parameter,'alpha')
                    obj.alpha = parameter.alpha;
                end
                
                if isfield(parameter,'te')
                    obj.te = parameter.te;
                end
                
                if isfield(parameter,'te_us')
                    obj.te_us = parameter.te_us;
                end
                
                obj.TE = obj.te + obj.te_us/1000;
                
                if isfield(parameter,'tr')
                    obj.tr = parameter.tr;
                end
                
                if isfield(parameter,'tr_extra_us')
                    obj.tr_extra_us = parameter.tr_extra_us;
                end
                
                obj.TR = obj.tr + obj.tr_extra_us/1000;
                
                if isfield(parameter,'pe1_order')
                    obj.pe1_order = parameter.pe1_order;
                end
                
                if isfield(parameter,'pe2_centric_on')
                    obj.pe2_centric_on = parameter.pe2_centric_on;
                end
                
                if isfield(parameter,'ti')
                    obj.ti = parameter.ti;
                end
                
                if isfield(parameter,'VFA_angles')
                    obj.VFA_angles = parameter.VFA_angles;
                end
                
                if isfield(parameter,'VFA_size')
                    obj.VFA_size = parameter.VFA_size;
                end
                
                if isfield(parameter,'frame_loop_on')
                    obj.frame_loop_on = parameter.frame_loop_on;
                end
                
                if isfield(parameter,'radial_on')
                    obj.radial_on = parameter.radial_on;
                end
                
                if isfield(parameter,'slice_nav')
                    obj.slice_nav = parameter.slice_nav;
                end
                
                if isfield(parameter,'no_samples_nav')
                    obj.no_samples_nav = parameter.no_samples_nav;
                end
                
                if isfield(parameter,'gp_var_mul')
                    obj.gp_var_mul = parameter.gp_var_mul;
                end
                
                if isfield(parameter,'date')
                    obj.date = parameter.date;
                end
                
                if isfield(parameter,'coil_scaling')
                    obj.coil_scaling = parameter.coil_scaling;
                end
                
                if isfield(parameter,'pixelshift1')
                    obj.pixelshift1 = parameter.pixelshift1;
                end
                
                if isfield(parameter,'pixelshift2')
                    obj.pixelshift2 = parameter.pixelshift2;
                end
                
                if isfield(parameter,'pe2_traj')
                    obj.pe2_traj = parameter.pe2_traj;
                end
                
                if isfield(parameter,'scanner')
                    obj.scanner = parameter.scanner;
                end
                
            end
            
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Check whether there are sufficient experiments to peform a valid reconstruction
        % ---------------------------------------------------------------------------------
        function [obj,message] = checkNumberOfExperiments(obj)
            
            if obj.EXPERIMENT_ARRAY * obj.NO_VIEWS * obj.NO_VIEWS_2 * obj.NO_SLICES < 1024
                obj.validData_flag = false;
                message = 'ERROR: Not enough k-lines for reconstruction ...';
            else
                message = '';
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Check the number of averages
        % ---------------------------------------------------------------------------------
        function [obj,message] = checkNumberOfAverages(obj)
            
            if obj.NO_AVERAGES > 1
                obj.validData_flag = false;
                message = 'ERROR: Number of averages should be 1 ...';
            else
                message = '';
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Check for multi-slab data
        % ---------------------------------------------------------------------------------
        function [obj,message] = checkForMultiSlab(obj)
            
            if strcmp(obj.dataType,'3D') && obj.NO_SLICES > 1
                obj.validData_flag = false;
                message = 'ERROR: Only 3D single-slab data supported ...';
            else
                message = '';
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Check for variable flip angle data
        % ---------------------------------------------------------------------------------
        function [obj,message] = checkForVFA(obj)
            
            if strcmp(obj.dataType,'3D') && obj.VFA_size > 1
                obj.vfaData_flag = true;
                obj = setVariableFlipAngles(obj);
                message = strcat('INFO:',{' '},num2str(obj.VFA_size),{' '},'flip angles detected ...');
            else
                message = '';
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Set the variable flip angles, sort the angles in groups
        % ---------------------------------------------------------------------------------
        function obj = setVariableFlipAngles(obj)
            
            a = obj.VFA_size;
            b = unique(obj.VFA_angles(1:obj.VFA_size),'Stable');
            c = length(b);
            d = obj.VFA_angles(1:obj.VFA_size);
            if a == c
                nr_flip_angles = a;     % FA1, FA2, FA3, FA4, .... = dym1, dyn2, dyn3, dyn4, ...
                ls_flip_angles = b;
            elseif mod(a,c) == 0
                nr_flip_angles = c;     % FA1, FA1, ..., FA2, FA2, ..., FA3, FA3, ... = dym1, dyn2, dyn3, dyn4, ...
                ls_flip_angles = b;
            else
                nr_flip_angles = a;     % each dynamic has its own flip-angle
                ls_flip_angles = d;
            end
            obj.VFA_size = nr_flip_angles;
            obj.VFA_angles = ls_flip_angles;
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Check if the data is a valid file with a motion navigator acquisition
        % ---------------------------------------------------------------------------------
        function obj = acquisitionType(obj)
            
            if obj.radial_on == 1
                
                obj.validData_flag = true;
                
                obj.dataType = '2Dradial';
                if ndims(obj.data{1}) == 3
                    for i=1:obj.nr_coils
                        obj.data{i} = permute(obj.data{i},[1,4,2,3]);
                    end
                end
                
                obj.primary_navigator_point = round(obj.NO_SAMPLES/2);
                obj.nr_repetitions = size(obj.data{1},1);
                
            elseif (obj.slice_nav == 1) && (obj.no_samples_nav > 0)
                
                obj.validData_flag = true;
                
                if obj.NO_VIEWS_2 > 1
                    
                    obj.dataType = '3D';                                        % 3D data
                    for i=1:obj.nr_coils
                        if obj.EXPERIMENT_ARRAY > 1
                            obj.data{i} = permute(obj.data{i},[1,3,2,4]);
                        else
                            obj.data{i} = permute(obj.data{i},[4,2,1,3]);
                        end
                    end
                    
                else
                    
                    obj.dataType = '2D';                                        % 2D single-slice data
                    if ndims(obj.data{1}) == 3
                        for i=1:obj.nr_coils
                            obj.data{i} = permute(obj.data{i},[1,4,2,3]);
                        end
                    end
                    
                    if obj.NO_SLICES > 1                                        % 2D multi-slice data
                        obj.dataType = '2Dms';
                        obj.multi2D_flag = true;
                    else
                        obj.multi2D_flag = false;
                    end
                    
                end
                
                obj.primary_navigator_point = obj.no_samples_nav;
                if obj.nr_nav_points_used > obj.primary_navigator_point
                    obj.nr_nav_points_used = obj.primary_navigator_point;
                end
                obj.nr_repetitions = size(obj.data{1},1);
                
            else
                
                obj.validData_flag = false;
                
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Guess which reconstruction type (systolic function, diastolic
        % function, scout, VFA
        % ---------------------------------------------------------------------------------
        function obj = guessRecoType(obj)
            
            if strcmp(obj.dataType,'2D') || strcmp(obj.dataType,'2Dms') || strcmp(obj.dataType,'2Dradial')
                
                if obj.nr_repetitions > 200
                    
                    obj.reco_guess = 'diastolic function';
                    
                else
                    
                    obj.reco_guess = 'systolic function';
                    
                end
                
                
            end
            
            if strcmp(obj.dataType,'2Dms')
                
                slices = [obj.s_angle_var ; obj.p_angle_var ; obj.r_angle_var]';
                nr_unique_slices = size(unique(slices,'rows'),1);
                
                if nr_unique_slices > 1
                    
                    obj.reco_guess = 'scout';
                    
                end
                
            end
            
            if strcmp(obj.dataType,'3D')
                
                obj.reco_guess = 'systolic function';
                
                if obj.vfaData_flag
                    
                    obj.reco_guess = 'variable flip-angle';
                    
                end
                
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Read the MRD footer
        % ---------------------------------------------------------------------------------
        function [obj,message] = readMRDfooter(obj,mrdfile)
            
            try
                
                % Read information from the header and footer first
                fid = fopen(mrdfile,'r');
                val = fread(fid,4,'int32');
                xdim = val(1);
                ydim = val(2);
                zdim = val(3);
                dim4 = val(4);
                fseek(fid,18,'bof');
                data_type=fread(fid,1, 'uint16');
                data_type = dec2hex(data_type);
                fseek(fid,152,'bof');
                val = fread(fid,2, 'int32');
                dim5 = val(1);
                dim6 = val(2);
                no_samples = xdim;
                no_views = ydim;
                no_views_2 = zdim;
                no_slices = dim4;
                no_echoes = dim5;
                no_expts = dim6;
                
                % Determine datatype
                if size(data_type,2)>1
                    onlydatatype = data_type(2);
                    iscomplex = 2;
                else
                    onlydatatype = data_type(1);
                    iscomplex = 1;
                end
                switch onlydatatype
                    case '0'
                        datasize = 1; % size in bytes
                    case '1'
                        datasize = 1; % size in bytes
                    case '2'
                        datasize = 2; % size in bytes
                    case '3'
                        datasize = 2; % size in bytes
                    case '4'
                        datasize = 4; % size in bytes
                    case '5'
                        datasize = 4; % size in bytes
                    case '6'
                        datasize = 8; % size in bytes
                    otherwise
                        datasize = 4; % size in bytes
                end
                
                % Fast forward to beginning of footer
                fseek(fid,512,'bof');
                num2read = no_expts*no_echoes*no_slices*no_views_2*no_views*no_samples*iscomplex;
                fseek(fid,num2read*datasize,'cof');
                
                % Read the footer
                obj.mrdFooter = char(fread(fid,Inf,'uchar')');
                fclose(fid);
                
                message = '';
                
            catch ME
                
                obj.validData_flag = false;
                message = ME.message;
                
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Read RPR file
        % ---------------------------------------------------------------------------------
        function [obj,message] = readRPRfile(obj,filename)
            
            try
                fid = fopen(filename,'r');
                obj.rprFile = char(fread(fid,Inf,'uchar')');
                fclose(fid);
                obj.rpr_flag = true;
                message = '';
            catch
                obj.rprFile = '';
                obj.rpr_flag = false;
                message = 'WARNING: rpr file not found ...';
            end
            
        end
        
        
        
        
        
    end
    
end

