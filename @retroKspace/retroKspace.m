classdef retroKspace
    
    % K-space data and sorting class for retrospective app
    
    properties
        
        % K-space data
        raw
        kSpace
        kSpaceMrd
        kSpaceAvg
        trajectory
        
        % Cardiac and respiratory binning
        cardBins
        respBins
        cardBinNrs
        respBinNrs
        
    end
    
    
    
    methods
        
        
        % ---------------------------------------------------------------------------------
        % Object constructor
        % ---------------------------------------------------------------------------------
        function obj = retroKspace
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Extract unsorted k-data from raw data
        % ---------------------------------------------------------------------------------
        function objKspace = extractData(objKspace, objData)
            
            objKspace.raw = cell(objData.nr_coils);
            
            switch objData.dataType
                
                case {'2D','2Dms','3D','3Dp2roud'}
                    for i = 1:objData.nr_coils
                        objKspace.raw{i} = objData.data{i}(:,:,:,objData.primary_navigator_point+objData.nr_nav_points_discarded+1:end);
                    end
                    
                case '2Dradial'
                    for i = 1:objData.nr_coils
                        objKspace.raw{i} = objData.data{i};
                    end
                    
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Assign bin times, relative binning
        % ---------------------------------------------------------------------------------
        function objKspace = assignBinTimes(objKspace, objNav, app)
            
            card_locations = objNav.heartTrigPoints;
            resp_locations = objNav.respTrigPoints;
            nr_card_frames = app.nrCardFrames;
            nr_resp_frames = app.nrRespFrames;
            
            nr_card = length(card_locations); % number of cardiac time points
            nr_resp = length(resp_locations); % number of respiratory time points
            
            % cardiac binning
            %
            %  INPUT: card_locations = array with the cardiac trigger points in units of samples
            %         nr_card_frames = number of desired cardiac bins
            %
            %  card_locations(i)                                card_locations(i+1)
            %       ||  j=1     j=2     j=3     j=4     j=5     j=nr_of_card_frames
            %       ||       |       |       |       |       |       ||
            %       ||       |       |       |       |       |       ||
            %       ||       |       |       |       |       |       ||
            %       cnt=1   cnt=2   cnt=3   cnt=4   cnt=5   cnt=6   cnt=7   cnt=...
            %      cbins(1) cbins(2) .....
            %
            %  RESULT: array with time-stamp of all the cardiac bins for all heartbeats in the measurement in units of samples
            %
            
            cbins = zeros(1,(nr_card-1)*nr_card_frames);
            
            cnt = 1;
            for i=1:nr_card-1
                for j=1:nr_card_frames
                    cbins(cnt)=card_locations(i)+(j-1)*(card_locations(i+1)-card_locations(i))/nr_card_frames;
                    cnt = cnt + 1;
                end
            end
            
            % respiratory binning
            
            rbins = zeros(1,(nr_resp-1)*nr_resp_frames);
            
            cnt = 1;
            for i=1:nr_resp-1
                for j=1:nr_resp_frames
                    rbins(cnt)=resp_locations(i)+(j-1)*(resp_locations(i+1)-resp_locations(i))/nr_resp_frames;
                    cnt = cnt + 1;
                end
            end
            
            objKspace.cardBins = cbins;
            objKspace.respBins = rbins;
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Assign bin times, absolute binning
        % ---------------------------------------------------------------------------------
        function objKspace = assignBinTimesAbs(objKspace, objNav, objData, app)
            
            card_locations = objNav.heartTrigPoints;
            resp_locations = objNav.respTrigPoints;
            heartrate = objNav.meanHeartRate;
            TR = objData.TR;
            nr_card_frames = app.nrCardFrames;
            nr_resp_frames = app.nrRespFrames;
            
            nr_card = length(card_locations); % number of cardiac time points
            nr_resp = length(resp_locations); % number of respiratory time points
            
            %  cardiac binning = ABSOLUTE BINNING
            %
            %  INPUT: card_locations = array with the cardiac trigger points in units of samples
            %         nr_card_frames = number of desired cardiac bins
            %
            %  card_locations(i)                                card_locations(i+1)
            %       ||  j=1     j=2     j=3     j=4     j=5     j=nr_of_card_frames
            %       ||       |       |       |       |       |       ||
            %       ||       |       |       |       |       |       ||
            %       ||       |       |       |       |       |       ||
            %       cnt=1   cnt=2   cnt=3   cnt=4   cnt=5   cnt=6   cnt=7   cnt=...
            %      cbins(1) cbins(2) .....
            %
            %  RESULT: cbins = array with time-stamp of all the cardiac bins for all heartbeats in the measurement in units of samples
            %
            
            % frame duration in sample tijd
            framedur = (60/heartrate)*(1000/TR)/nr_card_frames;
            
            cbins = zeros(2,(nr_card-1)*nr_card_frames);
            
            cnt = 1;
            for i=1:nr_card-1
                
                nrframes = floor((card_locations(i+1)-card_locations(i))/framedur); % number of frames for this particular heartbeat
                if nrframes>nr_card_frames
                    nrframes = nr_card_frames;
                end
                
                for j=1:nrframes
                    cbins(1,cnt) = card_locations(i)+(j-1)*framedur; % divide heart-beat in nrframes equal bins of duration framedur
                    cbins(2,cnt) = j;                                % immediately assign frame number
                    cnt = cnt + 1;
                end
                
            end
            
            cbins = cbins(:,1:cnt-1);
            
            % respiratory binning = regular relative binning
            
            rbins = zeros(1,(nr_resp-1)*nr_resp_frames);
            
            cnt = 1;
            for i=1:nr_resp-1
                for j=1:nr_resp_frames
                    rbins(cnt)=resp_locations(i)+(j-1)*(resp_locations(i+1)-resp_locations(i))/nr_resp_frames;
                    cnt = cnt + 1;
                end
            end
            
            objKspace.cardBins = cbins;
            objKspace.respBins = rbins;
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Assign bin frames, relative binning
        % ---------------------------------------------------------------------------------
        function objKspace = assignBinFrames(objKspace, objNav, objData, app)
            
            % assigns all the measured k-lines to a specific cardiac phase and respiratory phase bin
            
            bin_times_card = objKspace.cardBins';
            bin_times_resp = objKspace.respBins';
            resp_window = objNav.respWindow;
            nr_klines = objData.nrKlines;
            nr_card_frames = app.nrCardFrames;
            nr_resp_frames = app.nrRespFrames;
            
            startpoint = round(bin_times_resp(2)+1);                        % first measurement to be considered
            endpoint = round(bin_times_resp(length(bin_times_resp))-1);     % last measurement to be considered
            loc_maxj_card = length(bin_times_card);                         % last bin border
            loc_maxj_resp = length(bin_times_resp);
            
            % Assignmnent of all k-line acquisitions to a specific bin
            %
            % INPUT : bin_times_card = array of all time-stamps (units of samples) of the whole acquisition
            %
            %                       >    <--- decreases ------ bin_times_card(j)
            %                      ||            |
            %                      ||  j-th bin  |
            %                      ||            |
            %                      ||            |
            %       - increases -----> i-th time point / sample
            %
            % OUTPUT: For each time point / k-line sample assignment to cardiac bin (saw-tooth pattern)
            %
            %         /|      /|    /|        /|
            %        / |     / |   / |      /  |
            %       /  |   /   |  /  |    /    |  /     = PHASE BINNING, each heart-beat has equal number of bins
            %      /   |  /    | /   |  /      | /                       despite differences in heart beat duration
            %     /    |/      |/    |/        |/
            %     12345 1 23 45 12345 1 2 3 4 5  <- bins
            %
            
            card_assignments = zeros(nr_klines,1);                      % zero = no assignment (breathing, begin/end of data)
            parfor i=startpoint:endpoint                                % start search to which heartbeat the measurement belongs
  
                j=loc_maxj_card;
                while i<bin_times_card(j) || j==1 %#ok<*PFBNS> 
                    j=j-1;
                end
                
                card_assignments(i) = mod(j-1,nr_card_frames)+1;        % assign to bin frame number = j modulus nr_frames
                if nr_resp_frames==1 && resp_window(i)==1               % if measurement is during respiration and only 1 resp state, put back to 0 to discard this k-line
                        card_assignments(i) = 0;
                end
            end
            
            resp_assignments = zeros(nr_klines,1);                      % zero = no assignment (breathing, begin/end of data)
            parfor i=startpoint:endpoint                                % start search to which respiration the measurement belongs
                j=loc_maxj_resp;
                while i<bin_times_resp(j) || j==1
                    j=j-1;
                end
                resp_assignments(i) = mod(j-1,nr_resp_frames)+1;        % assign to bin frame number = j modulus nr_frames
            end
            
            objKspace.cardBinNrs = card_assignments;
            objKspace.respBinNrs = resp_assignments;
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Assign bin frames, absolute binning
        % ---------------------------------------------------------------------------------
        function objKspace = assignBinFramesAbs(objKspace, objNav, objData, app)
            
            bin_times_card = objKspace.cardBins';
            bin_times_resp = objKspace.respBins';
            resp_window = objNav.respWindow;
            nr_klines = objData.nrKlines;
            nr_resp_frames = app.nrRespFrames;
            
            % assigns all the measured k-lines to a specific cardiac phase and respiratory phase bin
            %
            % ABSOLUTE BINNING
            
            startpoint = round(bin_times_resp(2)+1);                        % first measurement to be considered
            endpoint = round(bin_times_resp(length(bin_times_resp))-1);     % last measurement to be considered
            loc_maxj_card = length(bin_times_card);                         % last bin border
            loc_maxj_resp = length(bin_times_resp);
            
            % Assignmnent of all k-line acquisitions to a specific bin
            %
            % INPUT : bin_times_card = array of all time-stamps (units of samples) of the whole acquisition
            %
            %                       >    <--- decreases ------ bin_times_card(j)
            %                      ||            |
            %                      ||  j-th bin  |
            %                      ||            |
            %                      ||            |
            %       - increases -----> i-th time point / sample
            %
            % OUTPUT: For each time point / k-line sample assignment to cardiac bin (saw-tooth pattern)
            %
            %
            %
            %         /|      /|    /|        /|
            %        / |     / |   / |      /  |
            %       /  |   /   |  /  |    /    |  /     = ABSOLUTE BINNING, time between frames is equal
            %      /   |  /    | /   |  /      | /                          despite differences in heart beat duration
            %     /    |/      |/    |/        |/                           -> later frames receive less data
            %     12345 1234567 12345 123456789  <- frames/bins
            %
           
            card_assignments = zeros(nr_klines,1);                          % zero = no assignment (breathing, begin/end of data)
            parfor i=startpoint:endpoint                                    % start search to which heartbeat the measurement belongs
                j=loc_maxj_card;
                
                while i < bin_times_card(j,1) || j==1
                    j=j-1;
                end
                
                card_assignments(i) = bin_times_card(j,2);                  % assign to bin frame number
                
                if nr_resp_frames==1 && resp_window(i)==1
                    card_assignments(i) = 0;                                % if measurement is during respiration and only 1 resp state, put back to 0 to discard this k-line
                end
                
            end
            
            % Respiratory binning is still done using relative binning
            
            resp_assignments = zeros(nr_klines,1);                          % zero = no assignment (breathing, begin/end of data)
            parfor i=startpoint:endpoint                                    % start search to which respiration the measurement belongs
                j=loc_maxj_resp;
                while i<bin_times_resp(j) || j==1
                    j=j-1;
                end
                resp_assignments(i) = mod(j-1,nr_resp_frames)+1;            % assign to bin frame number = j modulus nr_frames
            end
            
            objKspace.cardBinNrs = card_assignments;
            objKspace.respBinNrs = resp_assignments;
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Fill K-space 2D
        % ---------------------------------------------------------------------------------
        function objKspace = fillKspace2D(objKspace, objData, app)
            
            share = app.SharingEditField.Value;
            objKspace.kSpace = cell(objData.nr_coils);
            objKspace.kSpaceAvg = [];
            includeWindow = objData.includeWindow.*objData.excludeWindow;
            card_bin_ass = objKspace.cardBinNrs;
            resp_bin_ass = objKspace.respBinNrs;
            dimx = objData.dimx;
            dimy = objData.dimy;
            dimz = objData.dimz;
            nrKsteps = objData.nrKsteps;
            traj = objKspace.trajectory;
            nr_card_frames = app.nrCardFrames;
            nr_resp_frames = app.nrRespFrames;
            nrdynamics = app.nrDynamics;
            
            for coilnr = 1:objData.nr_coils
                
                rawdata = objKspace.raw{coilnr};
                
                % This function creates 3 arrays
                % (1) the kspace data sorted into the correct cardiac frames and phase-encoding positions
                % (2) an array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes
                % (3) an average cardiac navigator signal
                
                % Required input:
                %
                % rawdata               = unsorted k-space data
                % navheart              = navigator signal, used to construct an average heartbeat
                % nr_card_frames        = number of desired cardiac frames, must be consisted with bin assignments (bin_ass)
                % nr_resp_frames        = number of desired respiratory frames
                % nr_dynamics           = number of desired dynamics
                % dimz                  = number of slices
                % dimy                  = dimensions of the images: dimy (phase encoding)
                % nrKsteps              = number of k-lines in 1 repetition
                % dimx                  = dimensions of the images: dimx (readout)
                % card_bin_ass          = the cardiac bin assignment array for all measured k-lines
                % resp_bin_ass          = the respiratory bin assignment array for all measured k-lines
                % trajectory            = the k-space trajectory
                % includewindow         = data which should be include: 1 = yes, 0 = no
                % share                 = number of k-space sharing points
                
                sorted_kspace = complex(zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics));   % fill temp k-space with zeros
                sorted_averages = zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics);          % fill temp nr averages array with zeros
                nr_reps = size(rawdata,1);                                                                      % number of k-space repetitions
                unsorted_kspace = reshape(rawdata,[1,size(rawdata),1]);
                
                % Dynamics assignment
                totalk = nr_reps * nrKsteps * dimz;
                dyn_bin_ass = round(linspace(0.5, nrdynamics+0.49, totalk));       % list of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time
                
                % Sorting
                cnt = 0;
                
                for slice=1:dimz                    % loop over slices
                    
                    for i=1:nr_reps                 % loop through all repetitions
                        
                        for j=1:nrKsteps            % loop through all the phase-encoding steps
                            
                            cnt = cnt + 1;
                            
                            if (card_bin_ass(cnt) > 0) && (includeWindow(cnt) == 1)         % if assigment = 0, this acquisition is discarded
                                
                                kline = traj(mod(cnt - 1,nrKsteps) + 1);             % the phase-encoding step using the trajectory info
                                sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) = sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) + unsorted_kspace(1,i,slice,j,:,1);   % add the data to the correct k-position
                                sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) = sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) + 1;        % increase the number of averages with 1
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                % Temp k-space
                new_kspace = sorted_kspace;
                new_averages = sorted_averages;
                
                % Find center of k-space
                kspacesum = squeeze(sum(sorted_kspace,[1 2 3 6]));                 % sum over all slices frames and dynamics
                [row, col] = find(ismember(kspacesum, max(kspacesum(:))));         % coordinate of k-space maximum = center of k-space
                
                % Weighted view sharing
                if (share > 0) && (nr_card_frames > 1 || nr_resp_frames > 1)
                    
                    % respiratory or cardiac frames
                    nrframes = nr_card_frames;
                    if nr_resp_frames > 1
                        nrframes = nr_resp_frames;
                        new_kspace = permute(new_kspace,[2,1,3,4,5,6]);
                        new_averages = permute(new_averages,[2,1,3,4,5,6]);
                        sorted_kspace = permute(sorted_kspace,[2,1,3,4,5,6]);
                        sorted_averages = permute(sorted_averages,[2,1,3,4,5,6]);
                    end
                    
                    % determine share range
                    maxshare = round(max([nr_card_frames nr_resp_frames])/2);     % maximum number of shares
                    share(share > maxshare) = maxshare;
                    
                    % define ellipsoid regions
                    Ry = round(dimy/share/2);
                    Rx = round(dimx/share/2);
                    [Y,X] = ndgrid(1:dimy,1:dimx);
                    L = zeros(share,dimy,dimx);
                    for i = 1:share
                        L(i,:,:) = sqrt( ((row-Y)/(Ry*i)).^2 + ((col-X)/(Rx*i)).^2 ) <= 1;
                    end
                    C(1,:,:) = L(1,:,:);
                    if share > 1
                        for i = 2:share
                            C(i,:,:) = L(i,:,:) - L(i-1,:,:);
                        end
                    end
                    
                    % weights
                    weights = zeros(share);
                    for i = 1:share
                        for j = 1:share
                            weights(i,j) = retroKspace.gauss(i+j-1,share,0);
                        end
                    end
                    weights = 0.5*weights/max(weights(:));
                    
                    % apply sharing to k-space
                    for frame = 1:nrframes
                        
                        for i = -share:share
                            
                            sharedframe = frame + i;
                            sharedframe(sharedframe < 1) = nrframes - sharedframe - 1;
                            sharedframe(sharedframe > nrframes) = sharedframe - nrframes;
                            
                            if i~=0
                                
                                for j = 1:share
                                    
                                    ROI = reshape(squeeze(C(j,:,:)),[1 1 1 dimy dimx 1])*weights(j,abs(i));
                                    new_kspace(:,frame,:,:,:,:)   = new_kspace(:,frame,:,:,:,:)   + sorted_kspace(:,sharedframe,:,:,:,:)   .* ROI;
                                    new_averages(:,frame,:,:,:,:) = new_averages(:,frame,:,:,:,:) + sorted_averages(:,sharedframe,:,:,:,:) .* ROI;
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                    % respiratory or cardiac frames
                    if nr_resp_frames > 1
                        new_kspace = permute(new_kspace,[2,1,3,4,5,6]);
                        new_averages = permute(new_averages,[2,1,3,4,5,6]);
                    end
                    
                end
                
                % Normalize by number of averages
                new_kspace = new_kspace./new_averages;
                new_kspace(isnan(new_kspace)) = complex(0); % correct for NaN or Inf because of division by zero in case of missing k-lines
                new_kspace(isinf(new_kspace)) = complex(0);
                
                % Apply a circular Tukey filter
                filterwidth = 0.1;
                tukeyfilter(1,1,1,:,:,1) = retroKspace.circtukey2D(dimy,dimx,row,col,filterwidth);
                new_kspace = new_kspace.*tukeyfilter;
                
                % Report back
                objKspace.kSpace{coilnr} = new_kspace;
                objKspace.kSpaceAvg = new_averages;
                
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Fill K-space 3D
        % ---------------------------------------------------------------------------------
        function objKspace = fillKspace3D(objKspace, objData, app)
            
            share = app.SharingEditField.Value;
            objKspace.kSpace = cell(objData.nr_coils);
            objKspace.kSpaceAvg = [];
            includeWindow = objData.includeWindow.*objData.excludeWindow;
            card_bin_ass = objKspace.cardBinNrs;
            resp_bin_ass = objKspace.respBinNrs;
            dimx = objData.dimx;
            dimy = objData.dimy;
            dimz = objData.dimz;
            nrKsteps = objData.nrKsteps;
            traj = objKspace.trajectory;
            nr_card_frames = app.nrCardFrames;
            nr_resp_frames = app.nrRespFrames;
            nrdynamics = app.nrDynamics;
            
            for coilnr = 1:objData.nr_coils
                
                rawdata = objKspace.raw{coilnr};
                
                % This function creates 2 arrays
                % (1) the 3D kspace data sorted into the correct cardiac frames and phase-encoding positions
                % (2) an array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes
                
                % Required input:
                %
                % raw                   = unsorted k-space data
                % repstart              = first k-space repetition that is considered, used for partial reconstruction of the data
                % navheart              = navigator signal, used to construct an average heartbeat
                % nr_of_card_frames     = number of desired cardiac frames, must be consisted with bin assignments (bin_ass)
                % dimz                  = 2nd phase encoding dimension
                % dimy                  = dimensions of the images: dimy (phase encoding)
                % nrKsteps              = number of k-lines in 1 repetition
                % dimx                  = dimensions of the images: dimx (readout)
                % card_bin_ass          = the cardiac bin assignment array for all measured k-lines
                % resp_bin_ass          = the respiratory bin assignment array for all measured k-lines
                % traj                  = the k-space trajectory
                % includewindow         = data which should be include: 1 = yes, 0 = no
                % share                 = number of weighted view sharing points
                
                sorted_kspace = complex(zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics));    % fill temp k-space with zeros
                sorted_averages = zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics);  % fill temp nr averages array with zeros
                nr_reps = size(rawdata,1);                                                         % number of k-space repetitions
                unsorted_kspace = reshape(rawdata,[1,size(rawdata),1]);
                
                % dynamics assignment
                totalk = nr_reps * nrKsteps * dimz;
                dyn_bin_ass = round(linspace(0.5, nrdynamics+0.49, totalk));       % list of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time
                
                % adapt trajectory for 3D acqusition
                
                % the y-dimension
                traj3D_y = zeros(nr_reps * nrKsteps * dimz,1);
                cnt = 1;
                for i = 1:nr_reps
                    for j = 1:nrKsteps
                        for k = 1:dimz
                            traj3D_y(cnt) = traj(j);
                            cnt = cnt + 1;
                        end
                    end
                end
                
                % the z-dimension
                traj3D_z = zeros(nr_reps * nrKsteps * dimz,1);
                
                switch objData.pe2_centric_on
                    
                    case 0
                        
                        % linear in the 3rd dimension
                        cnt = 1;
                        for i = 1:nr_reps
                            for j = 1:nrKsteps
                                for k = 1:dimz
                                    traj3D_z(cnt) = k;  % linear
                                    cnt = cnt + 1;
                                end
                            end
                        end
                        
                    case 1
                        
                        % centric in the 3rd dimension
                        cnt = 1;
                        cf = retroKspace.centricFilling(dimz);
                        for i = 1:nr_reps
                            for j = 1:nrKsteps
                                for k = 1:dimz
                                    traj3D_z(cnt) = cf(k);  % centric
                                    cnt = cnt + 1;
                                end
                            end
                        end
                        
                    case 2
                        
                        % special case
                        cnt = 1;
                        for i = 1:nr_reps
                            for j = 1:nrKsteps
                                for k = 1:dimz
                                    traj3D_z(cnt) = objData.pe2_traj(k) + round(dimz/2) + 1;
                                    cnt = cnt + 1;
                                end
                            end
                        end
                        
                end
                
                % Do the filling of k-space
                cnt = 0;
                
                for i = 1:nr_reps                  % loop through all repetitions
                    
                    for j = 1:nrKsteps             % loop through all the phase-encoding steps
                        
                        for k = 1:dimz             % loop through phase-encoding 3rd dimension
                            
                            cnt = cnt + 1;
                            
                            if (card_bin_ass(cnt) > 0) && (includeWindow(cnt) == 1)     % if assigment == 0, this acquisition is discarded
                                
                                kline_y = traj3D_y(cnt);            % the phase-encoding step using the 3D trajectory info
                                kline_z = traj3D_z(cnt);            % the 2nd phase-encoding
                                
                                sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),kline_z,kline_y,:,dyn_bin_ass(cnt))   = sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),kline_z,kline_y,:,dyn_bin_ass(cnt))   + unsorted_kspace(1,i,k,j,:,1);       % add the data to the correct k-position
                                sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),kline_z,kline_y,:,dyn_bin_ass(cnt)) = sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),kline_z,kline_y,:,dyn_bin_ass(cnt)) + 1;                                  % increase the number of averages with 1
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                % temp new k-space copy
                new_kspace = sorted_kspace;
                new_averages = sorted_averages;
                
                % find center of k-space
                kspacesum = squeeze(sum(sorted_kspace,[1 2 6]));             % sum over all slices frames and dynamics
                [~,idx] = max(kspacesum(:));
                [lev, row, col] = ind2sub(size(kspacesum),idx);              % coordinate of k-space maximum = center of k-space
                
                % weighted view sharing
                if (share > 0) && (nr_card_frames > 1 || nr_resp_frames > 1)
                    
                    % respiratory of cardiac frames
                    nrframes = nr_card_frames;
                    if nr_resp_frames > 1
                        nrframes = nr_resp_frames;
                        new_kspace = permute(new_kspace,[2,1,3,4,5,6]);
                        new_averages = permute(new_averages,[2,1,3,4,5,6]);
                        sorted_kspace = permute(sorted_kspace,[2,1,3,4,5,6]);
                        sorted_averages = permute(sorted_averages,[2,1,3,4,5,6]);
                    end
                    
                    % determine share range
                    maxshare = round(max([nr_card_frames nr_resp_frames])/2); % maximum number of shares
                    share(share > maxshare) = maxshare;
                    weights = retroKspace.gauss(1:share+1,share,0);
                    weights = weights/max(weights);
                    
                    % define ellipsoid regions
                    Rz = round(dimz/share/2);
                    Ry = round(dimy/share/2);
                    Rx = round(dimx/share/2);
                    [Z,Y,X] = ndgrid(1:dimz,1:dimy,1:dimx);
                    L = zeros(share,dimz,dimy,dimx);
                    for i = 1:share
                        L(i,:,:,:) = sqrt( ((lev-Z)/(Rz*i)).^2 + ((row-Y)/(Ry*i)).^2 + ((col-X)/(Rx*i)).^2 ) <= 1;
                    end
                    C(1,:,:,:) = L(1,:,:,:);
                    if share > 1
                        for i = 2:share
                            C(i,:,:,:) = L(i,:,:,:) - L(i-1,:,:,:);
                        end
                    end
                    
                    % weights
                    for i = 1:share
                        for j = 1:share
                            weights(i,j) = retroKspace.gauss(i+j-1,share,0);
                        end
                    end
                    weights = 0.5*weights/max(weights(:));
                    
                    % apply sharing to k-space
                    for frame = 1:nrframes
                        
                        for i = -share:share
                            
                            sharedframe = frame + i;
                            sharedframe(sharedframe < 1) = nrframes - sharedframe - 1;
                            sharedframe(sharedframe > nrframes) = sharedframe - nrframes;
                            
                            if i~=0
                                
                                for j = 1:share
                                    
                                    ROI = reshape(squeeze(C(j,:,:,:)),[1 1 dimz dimy dimx 1])*weights(j,abs(i));
                                    new_kspace(:,frame,:,:,:,:)   = new_kspace(:,frame,:,:,:,:)   + sorted_kspace(:,sharedframe,:,:,:,:)   .* ROI;
                                    new_averages(:,frame,:,:,:,:) = new_averages(:,frame,:,:,:,:) + sorted_averages(:,sharedframe,:,:,:,:) .* ROI;
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                    % respiratory of cardiac frames
                    if nr_resp_frames > 1
                        new_kspace = permute(new_kspace,[2,1,3,4,5,6]);
                        new_averages = permute(new_averages,[2,1,3,4,5,6]);
                    end
                    
                end
                
                % normalize by number of averages
                new_kspace = new_kspace./new_averages;
                new_kspace(isnan(new_kspace)) = complex(0);     % correct for NaN because of division by zero in case of missing k-lines
                new_kspace(isinf(new_kspace)) = complex(0);
                
                % apply a circular Tukey filter
                filterwidth = 0.1;
                tukeyfilter(1,1,:,:,:,1) = retroKspace.circtukey3D(dimz,dimy,dimx,lev,row,col,filterwidth);
                new_kspace = new_kspace.*tukeyfilter;
                
                % report back
                objKspace.kSpace{coilnr} = new_kspace;
                objKspace.kSpaceAvg = new_averages;
                
            end
            
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Fill K-space 2D real-time
        % ---------------------------------------------------------------------------------
        function objKspace = fillKspace2D_RT(objKspace, objNav, objData, app)
            
            share = app.SharingEditField.Value;
            objKspace.kSpace = cell(objData.nr_coils);
            objKspace.kSpaceAvg = [];
            bin_times_card = objNav.heartTrigPoints;
            dimx = objData.dimx;
            dimy = objData.dimy;
            dimz = objData.dimz;
            nrKsteps = objData.nrKsteps;
            traj = objKspace.trajectory;
            nr_card_frames = app.nrCardFrames;
            nr_resp_frames = app.nrRespFrames;
           
            for coilnr = 1:objData.nr_coils
                
                rawdata = objKspace.raw{coilnr};
                
                % This function creates 3 arrays
                % (1) the kspace data sorted into the correct cardiac frames,  and phase-encoding positions
                % (2) an array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes
                % (3) the number of dynamics = number of heart-beats in the datasets
                
                % Required input:
                %
                % rawdata               = unsorted k-space data
                % navheart              = navigator signal, used to construct an average heartbeat
                % nr_of_card_frames     = number of desired cardiac frames, must be consisted with bin assignments (bin_ass)
                % nr_of_resp_frames     = number of desired respiratory frames
                % dimz                  = number of slices
                % dimy                  = dimensions of the images: dimy (phase encoding)
                % nrKsteps              = number of k-lines in 1 repetition
                % dimx                  = dimensions of the images: dimx (readout)
                % card_bin_ass          = the cardiac bin assignment array for all measured k-lines
                % traj                  = the k-space trajectory
                % share                 = weighted view sharing of neighboring data
                
                nrdynamics = round(length(bin_times_card)/nr_card_frames);                                     % number of dynamics equals the number of heartbeats in the acquisition
                sorted_kspace = complex(zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics));       % fill temp k-space with zeros
                sorted_averages = zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics);              % fill temp averages array with zeros
                nr_reps = size(rawdata,1);                                                                          % number of k-space repetitions
                nr_klines = nr_reps * nrKsteps * dimz;                                                         % total number of k-lines
                
                % unsorted k-lines
                rawdata = permute(rawdata,[3,1,2,4]);
                unsorted_klines = reshape(rawdata,[1,1,1,nr_klines,dimx]);
                
                % fill k-space
                fcnt = 1;
                
                for i = 1 : nrdynamics
                    
                    for j = 1 : nr_card_frames
                        
                        for k = 1:nr_klines
                            
                            if fcnt < length(bin_times_card)
                                
                                if k > bin_times_card(1,fcnt) && k < bin_times_card(1,fcnt+1)
                                    
                                    % the phase-encoding step using the trajectory info
                                    kline = traj(mod(k - 1,nrKsteps) + 1);
                                    
                                    % fill k-line
                                    sorted_kspace(1,j,1,kline,:,i) = sorted_kspace(1,j,1,kline,:,i) + unsorted_klines(1,1,1,k,:);   % add the data to the correct k-position
                                    sorted_averages(1,j,1,kline,:,i) = sorted_averages(1,j,1,kline,:,i) + 1;
                                    
                                end
                                
                            end
                            
                        end
                        
                        fcnt = fcnt + 1;
                        
                    end
                    
                end
                
                % find center of k-space
                kspacesum = squeeze(sum(sorted_kspace,[1 2 3 6]));                 % sum over all slices frames and dynamics
                [row, col] = find(ismember(kspacesum, max(kspacesum(:))));         % coordinate of k-space maximum = center of k-space
                
                % Temp k-space
                new_kspace = sorted_kspace;
                new_averages = sorted_averages;
                
                % Weighted view sharing
                if (share > 0) && (nrdynamics > 1)
                    
                    % determine share range
                    maxshare = 20;                                          % maximum number of shares
                    share(share > maxshare) = maxshare;
                    
                    % define ellipsoid regions
                    Ry = round(dimy/share/2);
                    Rx = round(dimx/share/2);
                    [Y,X] = ndgrid(1:dimy,1:dimx);
                    L = zeros(share,dimy,dimx);
                    for i = 1:share
                        L(i,:,:) = sqrt( ((row-Y)/(Ry*i)).^2 + ((col-X)/(Rx*i)).^2 ) <= 1;
                    end
                    C(1,:,:) = L(1,:,:);
                    if share > 1
                        for i = 2:share
                            C(i,:,:) = L(i,:,:) - L(i-1,:,:);
                        end
                    end
                    
                    % weights
                    weights = zeros(share);
                    for i = 1:share
                        for j = 1:share
                            weights(i,j) = retroKspace.gauss(i+j-1,share,0);
                        end
                    end
                    weights = 0.5*weights/max(weights(:));
                    
                    % apply sharing to k-space
                    for frame = 1:nrdynamics
                        
                        for i = -share:share
                            
                            sharedframe = frame + i;
                            
                            if sharedframe > 0 && sharedframe < nrdynamics
                                
                                if i~=0
                                    
                                    for j = 1:share
                                        
                                        ROI = reshape(squeeze(C(j,:,:)),[1 1 1 dimy dimx 1])*weights(j,abs(i));
                                        new_kspace(:,:,:,:,:,frame)   = new_kspace(:,:,:,:,:,frame)   + sorted_kspace(:,:,:,:,:,sharedframe)   .* ROI;
                                        new_averages(:,:,:,:,:,frame) = new_averages(:,:,:,:,:,frame) + sorted_averages(:,:,:,:,:,sharedframe) .* ROI;
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                % Normalize by number of averages
                new_kspace = new_kspace./new_averages;
                new_kspace(isnan(new_kspace)) = complex(0);
                new_kspace(isinf(new_kspace)) = complex(0);
                
                % Apply a circular Tukey filter
                filterwidth = 0.1;
                tukeyfilter(1,1,1,:,:,1) = retroKspace.circtukey2D(dimy,dimx,row,col,filterwidth);
                new_kspace = new_kspace.*tukeyfilter;
                
                % Report back
                objKspace.kSpace{coilnr} = new_kspace;
                objKspace.kSpaceAvg = new_averages;
                app.nrDynamics = nrdynamics;
                
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Fill K-space 3D P2ROUD
        % ---------------------------------------------------------------------------------
        function objKspace = fillKspace3D_P2ROUD(objKspace, objData, app)
            
            objKspace.kSpace = cell(objData.nr_coils);
            objKspace.kSpaceAvg = [];
            includeWindow = objData.includeWindow.*objData.excludeWindow;
            card_bin_ass = objKspace.cardBinNrs;
            resp_bin_ass = objKspace.respBinNrs;
            dimx = objData.dimx;
            dimy = objData.dimy;
            dimz = objData.dimz;
            nrKsteps = objData.nrKsteps;
            traj = objKspace.trajectory;
            nr_card_frames = app.nrCardFrames;
            nr_resp_frames = app.nrRespFrames;
            nrdynamics = app.nrDynamics;
            
            for coilnr = 1:objData.nr_coils
                
                rawdata = objKspace.raw{coilnr};
                
                % This function creates 2 arrays
                % (1) the 3D kspace data sorted into the correct cardiac frames and phase-encoding positions
                % (2) an array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes
                % (3) an average cardiac navigator signal
                
                % Required input:
                %
                % raw                   = unsorted k-space data
                % repstart              = first k-space repetition that is considered, used for partial reconstruction of the data
                % navheart              = navigator signal, used to construct an average heartbeat
                % nr_of_card_frames     = number of desired cardiac frames, must be consisted with bin assignments (bin_ass)
                % dimz                  = 2nd phase encoding dimension
                % dimy                  = dimensions of the images: dimy (phase encoding)
                % nrKsteps              = number of k-lines in 1 repetition
                % dimx                  = dimensions of the images: dimx (readout)
                % card_bin_ass          = the cardiac bin assignment array for all measured k-lines
                % resp_bin_ass          = the respiratory bin assignment array for all measured k-lines
                % traj                  = the k-space trajectory
                % includewindow         = data which should be include: 1 = yes, 0 = no
                
                sorted_kspace = complex(zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics));       % fill temp k-space with zeros
                sorted_averages = zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics);              % fill temp nr averages array with zeros
                nr_reps = size(rawdata,1);                                                                      % number of k-space repetitions
                unsorted_kspace = reshape(rawdata,[1,size(rawdata),1]);
                
                % dynamics assignment
                totalk = nr_reps * nrKsteps * dimz;
                dyn_bin_ass = round(linspace(0.5, nrdynamics+0.49, totalk));       % list of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time
                
                % trajectory for 3D p2roud acquisition
                offset1 = floor(dimy/2 + 1);
                offset2 = floor(dimz/2 + 1);
                cnt1 = 1;
                for i = 1:nr_reps
                    cnt2 = 1;
                    for j = 1:nrKsteps
                        for k = 1:dimz
                            traj3D_y(cnt1) = traj(cnt2) + offset1; %#ok<*AGROW> 
                            traj3D_z(cnt1) = traj(cnt2+1) + offset2;
                            cnt1 = cnt1 + 1;
                            cnt2 = cnt2 + 2;
                        end
                    end
                end
                
                % do the filling of k-space
                cnt = 0;
                
                for i = 1:nr_reps                  % loop through all repetitions
                    
                    for j = 1:nrKsteps             % loop through all the phase-encoding steps
                        
                        for k = 1:dimz             % loop through phase-encoding 3rd dimension
                            
                            cnt = cnt + 1;
                            
                            if (card_bin_ass(cnt) > 0) && (includeWindow(cnt) == 1)     % if assigment == 0, this acquisition is discarded
                                
                                kline_y = traj3D_y(cnt);            % the phase-encoding step using the 3D trajectory info
                                kline_z = traj3D_z(cnt);            % the 2nd phase-encoding
                                
                                sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),kline_z,kline_y,:,dyn_bin_ass(cnt)) = sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),kline_z,kline_y,:,dyn_bin_ass(cnt)) + unsorted_kspace(1,i,k,j,:,1);     % add the data to the correct k-position
                                sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),kline_z,kline_y,:,dyn_bin_ass(cnt)) = sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),kline_z,kline_y,:,dyn_bin_ass(cnt)) + 1;                            % increase the number of averages with 1
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                % normalize by number of averages
                sorted_kspace = sorted_kspace./sorted_averages;       % normalize by number of averages
                sorted_kspace(isnan(sorted_kspace)) = complex(0);     % correct for NaN because of division by zero in case of missing k-lines
                sorted_kspace(isinf(sorted_kspace)) = complex(0);
                
                % apply a circular Tukey filter
                filterwidth = 0.1;
                tukeyfilter(1,1,:,:,:,1) = retroKspace.circtukey3D(dimz,dimy,dimx,lev,row,col,filterwidth);
                sorted_kspace = sorted_kspace.*tukeyfilter;
                
                % report back
                objKspace.kSpace{coilnr} = sorted_kspace;
                objKspace.kSpaceAvg = sorted_averages;
                
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Fill K-space RADIAL
        % ---------------------------------------------------------------------------------
        function objKspace = fillKspace2D_RADIAL(objKspace, objNav, objData)
            
            objKspace.kSpace = cell(objData.nr_coils);
            objKspace.kSpaceAvg = [];
            includeWindow = objData.includeWindow.*objData.excludeWindow;
            card_bin_ass = objKspace.cardBinNrs;
            resp_bin_ass = objKspace.respBinNrs;
            dimx = objData.dimx;
            dimy = objData.dimy;
            dimz = objData.dimz;
            nrKsteps = objData.nrKsteps;
            traj = objKspace.trajectory;
            nr_card_frames = app.nrCardFrames;
            nr_resp_frames = app.nrRespFrames;
            nrdynamics = app.nrDynamics;
            
            for coilnr = 1:objData.nr_coils
                
                rawdata = objKspace.raw{coilnr};
                phaseoffset = objNav.nav_phase{coilnr};
                
                % This function creates 3 arrays
                % (1) the kspace data sorted into the correct cardiac frames and phase-encoding positions
                % (2) an array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes
                % (3) an average cardiac navigator signal
                
                % Required input:
                %
                % rawdata               = unsorted k-space data
                % navheart              = navigator signal, used to construct an average heartbeat
                % nr_card_frames        = number of desired cardiac frames, must be consisted with bin assignments (bin_ass)
                % nr_resp_frames        = number of desired respiratory frames
                % nr_of_dynamics        = number of desired dynamics
                % dimz                  = number of slices
                % dimy                  = dimensions of the images: dimy (phase encoding)
                % nrKsteps              = number of k-lines in 1 repetition
                % dimx                  = dimensions of the images: dimx (readout)
                % card_bin_ass          = the cardiac bin assignment array for all measured k-lines
                % resp_bin_ass          = the respiratory bin assignment array for all measured k-lines
                % traj                  = the k-space trajectory
                % includewindow         = data which should be include: 1 = yes, 0 = no
                
                sorted_kspace = complex(zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics));   % fill temp k-space with zeros
                sorted_averages = zeros(nr_resp_frames,nr_card_frames,dimz,dimy,dimx,nrdynamics);          % fill temp nr averages array with zeros
                nr_reps = size(rawdata,1);                                                                 % number of k-space repetitions
                unsorted_kspace = reshape(rawdata,[1,size(rawdata),1]);
                
                % dynamics assignment
                totalk = nr_reps * nrKsteps * dimz;
                dyn_bin_ass = round(linspace(0.5, nrdynamics+0.49, totalk));       % list of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time
                
                cnt = 0;
                
                for slice=1:dimz                    % loop over slices
                    
                    for r=1:nr_reps                 % loop through all repetitions
                        
                        for y=1:nrKsteps           % loop through all the phase-encoding steps
                            
                            cnt = cnt + 1;
                            
                            if (card_bin_ass(cnt) > 0) && (includeWindow(cnt) == 1)     % if assigment = 0, this acquisition is discarded
                                
                                kline = traj(mod(cnt - 1,nrKsteps) + 1);         % the phase-encoding step using the trajectory info
                                sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) = sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) + unsorted_kspace(1,r,slice,y,:,1)*exp(-1i*phaseoffset(cnt));   % add the data to the correct k-position
                                sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) = sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) + 1;        % increase the number of averages with 1
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                % find center of k-space
                kspacesum = squeeze(sum(sorted_kspace,[1 2 3 6]));                  % sum over all slices frames and dynamics
                [row, col] = find(ismember(kspacesum, max(kspacesum(:))));          % coordinate of k-space maximum = center of k-space
                
                sorted_kspace = sorted_kspace./sorted_averages;                     % normalize by number of averages
                sorted_kspace(isnan(sorted_kspace)) = complex(0);                   % correct for NaN or Inf because of division by zero in case of missing k-lines
                sorted_kspace(isinf(sorted_kspace)) = complex(0);
          
                % imshow(squeeze(real(sorted_kspace(1,1,1,:,:))),[]);
                % plot(squeeze(abs(kspace(1,1,1,:,64))));
                
                % Apply a circular Tukey filter
                filterwidth = 0.1;
                tukeyfilter(1,1,1,:,:,1) = retroKspace.circtukey2D(dimy,dimx,row,col,filterwidth);
                sorted_kspace = sorted_kspace.*tukeyfilter;
                
                % Report back
                objKspace.kSpace{coilnr} = sorted_kspace;
                objKspace.kSpaceAvg = sorted_averages;
                
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Reshape k-space
        % ---------------------------------------------------------------------------------
        function objKspace = reshapeKspace(objKspace)
            
            % Reshape to either cardiac or respiratory CINE
            
            [s1,s2,s3,s4,s5,s6] = size(objKspace.kSpaceAvg);
            s = max([s1,s2]);
            nr_coils = length(objKspace.kSpace);
            for i = 1:nr_coils
                objKspace.kSpace{i} = reshape(objKspace.kSpace{i},[s,s3,s4,s5,s6]);
                objKspace.kSpace{i} = permute(objKspace.kSpace{i},[1,3,4,2,5,6]);
            end
            objKspace.kSpaceAvg = reshape(objKspace.kSpaceAvg,[s,s3,s4,s5,s6]);
            objKspace.kSpaceAvg = permute(objKspace.kSpaceAvg,[1,3,4,2,5,6]);
            
        end
        

        
        % ---------------------------------------------------------------------------------
        % K-space statistics
        % ---------------------------------------------------------------------------------
        function [obj, app] = kSpaceStats(obj, app)
            
            app.FillingViewField.Value = round(100*nnz(obj.kSpaceAvg)/numel(obj.kSpaceAvg));
            if app.FillingViewField.Value < 20
                app.SetStatus(1);
                app.TextMessage('WARNING: Low k-space filling, decrease number of frames or dynamics ...');
            end

            % Warning in case of large dataset
            if app.FramesEditField.Value*app.DynamicsEditField.Value > 2000
                app.SetStatus(1);
                app.TextMessage('WARNING: large dataset, reconstruction may take a very long time ...');
            end

        end % kSpaceStats


        
        % ---------------------------------------------------------------------------------
        % K-space trajectory
        % ---------------------------------------------------------------------------------
        function [objKspace, objData] = getKspaceTrajectory(objKspace, objData, app)
            
            % Trajectory (linear, zigzag, user-defined, radial)
            switch objData.pe1_order
                
                case 0
                    objKspace.trajectory = linear_trajectory(objData.nrKsteps);
                    
                case 2
                    objKspace.trajectory = zigzag_trajectory(objData.nrKsteps);
                    
                case 3
                    objKspace.trajectory = array_trajectory(objData.NO_VIEWS,objData.gp_var_mul);
                    
                case 4 % reserved for 3D P2ROUD acquisition
                    flist = dir(fullfile(app.mrd_import_path,'*.txt'));
                    if ~isempty(flist)
                        objKspace.trajectory = load([flist(1).folder,filesep,flist(1).name]);
                        objData.dataType = '3Dp2roud';
                        app.TextMessage(strcat('P2ROUD trajectory: ',{' '},flist(1).name));
                    else
                        objData.dataType = '3D';
                        objKspace.trajectory = linear_trajectory(objData.nrKsteps);
                        app.TextMessage('WARNING: P2ROUD trajectory file not found, assuming linear k-space filling ...');
                        app.SetStatus(1);
                    end
                    
                case 8 % reserved for radial
                    objKspace.trajectory = radial_trajectory(objData.nrKsteps);
                    
            end
            
            % ---------------------------------------------------------------------------------
            % Linear trajectory
            % ---------------------------------------------------------------------------------
            function traject = linear_trajectory(nr_pe)
                
                for i=1:nr_pe
                    traj(i)=i;
                end
                
                traject=traj;
                
            end
            
            % ---------------------------------------------------------------------------------
            % Zig-zag trajectory
            % ---------------------------------------------------------------------------------
            function traject = zigzag_trajectory(nr_pe)
                
                for i=1:nr_pe/2
                    traj(i)=2*(i-1)+1;
                    traj(i+nr_pe/2)=nr_pe-2*i+2;
                end
                
                traject=traj;
                
            end
            
            % ---------------------------------------------------------------------------------
            % Array trajectory
            % ---------------------------------------------------------------------------------
            function traject = array_trajectory(nr_pes,pe_array)
                
                traject = round(pe_array(1:nr_pes) - min(pe_array) + 1);
                
            end
            
            % ---------------------------------------------------------------------------------
            % Radial trajectory
            % ---------------------------------------------------------------------------------
            function traject = radial_trajectory(nr_pe)
                
                % Radial k-space trajectory
                % The actual trajectory is later incorporated in the trajectory that is used by the Bart toolbox
                % Therefore the trajectory is simply linear here
                
                for i=1:nr_pe
                    traject(i) = i;
                end
                
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Back to K-space
        % ---------------------------------------------------------------------------------
        function objKspace = backToKspace(objKspace, objData, objReco)

            objKspace.kSpaceMrd = [];
            [nf, ~, ~, dimz, nr] = size(objReco.movieExp);

            switch objData.dataType

                case {'2D','2Dms','radial'}

                    for i = 1:nf
                        for j = 1:nr
                            for k = 1:dimz
                                objKspace.kSpaceMrd(i,:,:,k,j) = retroKspace.fft2r(squeeze(objReco.movieExp(i,:,:,k,j)));
                            end
                        end
                    end

                    % samples, views, views2, slices, echoes (frames), experiments
                    objKspace.kSpaceMrd = permute(objKspace.kSpaceMrd,[2,3,6,4,1,5]);
    
                case {'3D','3Dp2roud'}

                    for i = 1:nf
                        for j = 1:nr
                            objKspace.kSpaceMrd(i,:,:,:,j) = retroKspace.fft3r(squeeze(objReco.movieExp(i,:,:,:,j)));
                        end
                    end

                    % samples, views, views2, slices, echoes (frames), experiments
                    objKspace.kSpaceMrd = permute(objKspace.kSpaceMrd,[2,3,4,6,1,5]);

            end

        end


    end % Public methods
    
    
    
    % ---------------------------------------------------------------------------------
    % Static methods
    % ---------------------------------------------------------------------------------
    methods (Static)
        
        
        % ---------------------------------------------------------------------------------
        % 2D Tukey filter
        % ---------------------------------------------------------------------------------
        function output = circtukey2D(dimy, dimx, row, col, filterwidth)
            
            domain = 256;
            base = zeros(domain,domain);
            
            tukey1 = tukeywin(domain,filterwidth);
            tukey1 = tukey1(domain/2+1:domain);
            
            shifty = (row-dimy/2)*domain/dimy;
            shiftx = (col-dimx/2)*domain/dimx;
            
            y = linspace(-domain/2, domain/2, domain);
            x = linspace(-domain/2, domain/2, domain);
            
            for i=1:domain
                
                for j=1:domain
                    
                    rad = round(sqrt((shiftx-x(i))^2 + (shifty-y(j))^2));
                    
                    if (rad <= domain/2) && (rad > 0)
                        
                        base(j,i) = tukey1(rad);
                        
                    end
                    
                end
                
            end
            
            output = imresize(base,[dimy dimx]);
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % 3D Tukey filter
        % ---------------------------------------------------------------------------------
        function output = circtukey3D(dimz,dimy,dimx,lev,row,col,filterwidth)
            
            domain = 256;
            
            base = zeros(domain,domain,domain);
            
            tukey1 = tukeywin(domain,filterwidth);
            tukey1 = tukey1(domain/2+1:domain);
            
            shiftz = (lev-dimz/2)*domain/dimz;
            shifty = (row-dimy/2)*domain/dimy;
            shiftx = (col-dimx/2)*domain/dimx;
            
            z = linspace(-domain/2, domain/2, domain);
            y = linspace(-domain/2, domain/2, domain);
            x = linspace(-domain/2, domain/2, domain);
            
            for i=1:domain
                
                for j=1:domain
                    
                    for k = 1:domain
                        
                        rad = round(sqrt((shiftx-x(i))^2 + (shifty-y(j))^2 + (shiftz-z(k))^2));
                        
                        if (rad <= domain/2) && (rad > 0)
                            
                            base(k,j,i) = tukey1(rad);
                            
                        end
                        
                    end
                    
                end
                
            end
            
            output = imresize3(base,[dimz dimy dimx]);
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Centric K-space filling scheme
        % ---------------------------------------------------------------------------------
        function scheme = centricFilling(no_views_2)
            
            ord2 = zeros(2*round(no_views_2/2));
            
            for g=1:round(no_views_2/2)
                
                ord2(2*g-1)=no_views_2/2+g;
                ord2(2*g)=no_views_2/2-g+1;
                
            end % centricFilling
            
            scheme = round(ord2);
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Gauss function
        % ---------------------------------------------------------------------------------
        function y = gauss(x,s,m)
            
            % GAUSS  Gaussian function
            %
            %  Y = GAUSS( X , S , M )
            %
            %  Y = EXP(-(X-M).^2./S.^2)./(sqrt(2*pi).*S);
            %
            %  sum( Y(X=(-inf..inf)) * dX ) = 1/sqrt(2)
            
            Nin = nargin;
            
            if Nin < 2
                s = 1;
            end
            if Nin < 3
                m = 0;
            end
            
            x = ((x-m).^2) ./ (s.^2);
            
            s = sqrt(2*pi) * s;
            
            y = exp(-x) ./ s;
            
        end
        


        % ---------------------------------------------------------------------------------
        % 2D FFT
        % ---------------------------------------------------------------------------------
        function Y = fft2r(x)

            Y = fftshift(ifft(fftshift(x,1),[],1),1)*sqrt(size(x,1));
            Y = fftshift(ifft(fftshift(Y,2),[],2),2)*sqrt(size(x,2));
                        
        end



        % ---------------------------------------------------------------------------------
        % 3D FFT
        % ---------------------------------------------------------------------------------
        function Y = fft3r(x)

            Y = fftshift(ifft(fftshift(x,1),[],1),1)*sqrt(size(x,1));
            Y = fftshift(ifft(fftshift(Y,2),[],2),2)*sqrt(size(x,2));
            Y = fftshift(ifft(fftshift(Y,3),[],3),3)*sqrt(size(x,3));

        end




    end % methods
    
end % retroKspace

