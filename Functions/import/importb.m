function [rawdata,parameters] = importb(import_path)


% Parameters
info1 = jcampread(strcat(import_path,'acqp'));
info2 = jcampread(strcat(import_path,'method'));


% Scanner type
parameters.scanner = 'B-type';


% Slices
parameters.NO_SLICES = str2num(info1.NSLICES);
parameters.SLICE_THICKNESS = str2num(info2.pvm.slicethick) * parameters.NO_SLICES;



% Matrix in readout direction
parameters.NO_SAMPLES = info1.acq.size(1) / 2;
if isfield(info2.pvm,"matrix")
    parameters.NO_VIEWS = info2.pvm.encmatrix(1);
end



% Matrix in phase encoding direction
parameters.NO_VIEWS = info1.acq.size(2);
if isfield(info2.pvm,"matrix")
    parameters.NO_VIEWS = info2.pvm.encmatrix(2);
end



% Phase encoding orientation
parameters.PHASE_ORIENTATION = 1;
pm1 = -1;
pm2 = -1;
if isfield(info2.pvm,'spackarrreadorient')
   if strcmp(info2.pvm.spackarrreadorient,'L_R')
      parameters.PHASE_ORIENTATION = 0;  
      flr =  1;
      pm1 = +1;
      pm2 = -1;
   end
   if strcmp(info2.pvm.spackarrreadorient,'A_P')
      parameters.PHASE_ORIENTATION = 1;  
      flr =  0;
      pm1 = -1;
      pm2 = -1;
   end
   if strcmp(info2.pvm.spackarrreadorient,'H_F')
      parameters.PHASE_ORIENTATION = 1;  
      flr =  0;
      pm1 = -1;
      pm2 = -1;
   end
end



% Matrix in 2nd phase encoding direction
parameters.NO_VIEWS_2 = 1;
parameters.pe2_centric_on = 0;



% FOV
parameters.FOV = info1.acq.fov(1)*10;
parameters.FOV2 = info1.acq.fov(2)*10;
parameters.FOVf = round(8*parameters.FOV2/parameters.FOV);



% Sequence parameters
parameters.tr = info1.acq.repetition_time;
parameters.te = info1.acq.echo_time;
parameters.alpha = str2num(info1.acq.flip_angle);
parameters.NO_ECHOES = 1;
parameters.NO_AVERAGES = str2num(info1.NA);



% Other parameters
parameters.date = datetime;
parameters.nucleus = '1H';
parameters.PPL = 'Navigator Sequence';
parameters.filename = 'Proton';
parameters.field_strength = str2num(info1.BF1)/42.58;
parameters.filename = 111;
parameters.pe1_order = 2;
parameters.radial_on = 0;
parameters.slice_nav = 1;


% Number of navigator points
if isfield(info2.pvm,"navpoints")
    parameters.no_samples_nav = str2num(info2.pvm.navpoints);
else
    parameters.no_samples_nav = str2num(info2.NavSize) / 2;
end



% Number of receiver coils
parameters.nr_coils = str2num(info2.pvm.encnreceivers);



% Trajectory 1st phase encoding direction
if isfield(info2.pvm,'ppggradamparray1')
    if isfield(info2.pvm,'enczfaccel1') && isfield(info2.pvm,'encpftaccel1')
        parameters.gp_var_mul = round(pm1 * info2.pvm.ppggradamparray1 * str2num(info2.pvm.enczfaccel1) * str2num(info2.pvm.encpftaccel1) * (parameters.NO_VIEWS / 2 - 0.5));
    else
        parameters.gp_var_mul = round(pm1 * info2.pvm.ppggradamparray1 * (parameters.NO_VIEWS / 2 - 0.5));
    end
    parameters.pe1_order = 3;
elseif isfield(info2.pvm,'encvalues1')
    if isfield(info2.pvm,'enczf') && isfield(info2.pvm,'encpft')
        parameters.gp_var_mul = round(pm1 * info2.pvm.encvalues1 * info2.pvm.enczf(2) * info2.pvm.encpft(2) * (parameters.NO_VIEWS / 2 - 0.5));
    else
        parameters.gp_var_mul = round(pm1 * info2.pvm.encvalues1 * (parameters.NO_VIEWS / 2 - 0.5));
    end
    parameters.pe1_order = 3;
else
    % assume zigzag
    parameters.pe1_order = 2;
end



% Data type
datatype = 'int32';
if isfield(info1.acq,'word_size')
    if strcmp(info1.acq.word_size,'_32_BIT') 
        datatype = 'int32';
    end
     if strcmp(info1.acq.word_size,'_16_BIT') 
        datatype = 'int16';
    end
end



% Read data
if isfile(strcat(import_path,'fid.orig'))
    fileID = fopen(strcat(import_path,'fid.orig'));
else
    fileID = fopen(strcat(import_path,'rawdata.job0'));
end
data = fread(fileID,datatype);
fclose(fileID);
kreal = data(1:2:end);
kim = data(2:2:end);
kspace = kreal + 1j*kim;



% Read navigator
if isfile(strcat(import_path,'fid.NavFid'))
    fileID = fopen(strcat(import_path,'fid.NavFid'));
else
    fileID = fopen(strcat(import_path,'rawdata.job1'));
end
navdata = fread(fileID,datatype);
fclose(fileID);
kreal = navdata(1:2:end);
kim = navdata(2:2:end);
navkspace = kreal + 1j*kim;



% Phase offset
if isfield(info1.acq,'phase1_offset')
    parameters.pixelshift1 = round(pm1 * parameters.NO_VIEWS * info1.acq.phase1_offset / parameters.FOV);
end




% -----------------------------------------------------------------------------------------------
%                   2D DATA
% -----------------------------------------------------------------------------------------------

if strcmp(info2.pvm.spatdimenum,"2D") || strcmp(info2.pvm.spatdimenum,"<2D>")

    % imaging k-space
    kspace = reshape(kspace,parameters.NO_SLICES,parameters.NO_SAMPLES,parameters.nr_coils,parameters.NO_VIEWS,[]);
    parameters.EXPERIMENT_ARRAY = size(kspace,5);
    kspace = permute(kspace,[3,5,1,4,2]); % nc, nr, ns, np, nf
    
    % flip readout if needed
    if flr
       kspace = flip(kspace,5); 
    end
   
    % coil intensity scaling
    if isfield(info2.pvm,'encchanscaling')
        for i = 1:parameters.nr_coils
            kspace(i,:) = kspace(i,:) * info2.pvm.encchanscaling(i);
        end
    end
    
    % navigator
    navkspace = reshape(navkspace,parameters.NO_SLICES,parameters.no_samples_nav,parameters.nr_coils,parameters.NO_VIEWS,parameters.EXPERIMENT_ARRAY);
    navkspace = permute(navkspace,[3,5,1,4,2]);
    
    % 34 point spacer
    kspacer = zeros(parameters.nr_coils,parameters.EXPERIMENT_ARRAY,parameters.NO_SLICES,parameters.NO_VIEWS,34);
    
    % combine navigator + spacer + k-space
    raw = cat(5,navkspace,kspacer,kspace);
    for i = 1:parameters.nr_coils
        rawdata{i} = squeeze(raw(i,:,:,:,:));
        rawdata{i} = reshape(rawdata{i},parameters.EXPERIMENT_ARRAY,parameters.NO_SLICES,parameters.NO_VIEWS,parameters.NO_SAMPLES+34+parameters.no_samples_nav);
    end
    
end




% -----------------------------------------------------------------------------------------------
%                   3D DATA (NOT REALLY TESTED)
% -----------------------------------------------------------------------------------------------

if strcmp(info2.pvm.spatdimenum,"3D") || strcmp(info2.pvm.spatdimenum,"<3D>")
    
    % 2nd phase encoding direction
    parameters.NO_VIEWS_2 = info1.acq.size(3);
    if isfield(info2.pvm,"matrix")
        parameters.NO_VIEWS = info2.pvm.encmatrix(3);
    end
    
    % phase offset 2
    if isfield(info1.acq,'phase2_offset')
        parameters.pixelshift2 = round(pm2 * parameters.NO_VIEWS_2 * info1.acq.phase2_offset/parameters.FOV);
    end
    
    % slice thickness
    parameters.SLICE_THICKNESS = str2num(info2.pvm.slicethick);
    
    % 2nd phase encoding trajectory
    parameters.pe2_centric_on = 0;
    if isfield(info2.pvm,"encsteps2")
        parameters.pe2_traj = info2.pvm.encsteps2;
        parameters.pe2_centric_on = 2;
    end
    if isfield(info2.pvm,'encvalues2')
        parameters.pe2_traj = round(info2.pvm.encvalues2 * (parameters.NO_VIEWS_2/2-0.5));
        parameters.pe2_centric_on = 2;
    end
    
    % k-space
    kspace = reshape(kspace,parameters.nr_coils,parameters.NO_SAMPLES,parameters.NO_VIEWS,parameters.NO_VIEWS_2,[]);
    parameters.EXPERIMENT_ARRAY = size(kspace,5);
    kspace = permute(kspace,[1,5,4,3,2]);
    
    % flip readout if needed
    if flr
       kspace = flip(kspace,5); 
    end
    
    % coil intesnity scaling
    if isfield(info2.pvm,'encchanscaling')
        for i = 1:parameters.nr_coils
            kspace(i,:) = kspace(i,:) * info2.pvm.encchanscaling(i);
        end
    end
    
    % navigator
    navkspace = reshape(navkspace,parameters.nr_coils,parameters.no_samples_nav,parameters.NO_VIEWS,parameters.NO_VIEWS_2,parameters.EXPERIMENT_ARRAY);
    navkspace = permute(navkspace,[1,5,4,3,2]);
    
    % 34 point spacer
    kspacer = zeros(parameters.nr_coils,parameters.EXPERIMENT_ARRAY,parameters.NO_VIEWS_2,parameters.NO_VIEWS,34);
    
    % combine navigator + spacer + k-space
    raw = cat(5,navkspace,kspacer,kspace);
    for i = 1:parameters.nr_coils
        rawdata{i} = squeeze(raw(i,:,:,:,:));
        rawdata{i} = reshape(rawdata{i},parameters.EXPERIMENT_ARRAY,parameters.NO_VIEWS_2,parameters.NO_VIEWS,parameters.NO_SAMPLES+34+parameters.no_samples_nav);
    end
    
end





