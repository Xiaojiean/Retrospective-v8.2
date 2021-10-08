classdef retroNav
    
    % Navigator data class for retrospective app
    
    properties
      
        % Raw navigator
        nav_amplitude
        nav_phase
        up_down = 1
        
        % Filtering
        physiofiltersettings
        detectedHR
        detectedRR
        powerspectrum
        frequency
        bandwidthHR
        bandwidthRR
        heartNav
        respNav
        PCANav
        
        % Triggering assignments
        splineFactor = 60           % data interpolation factor to prevent navigator discretization by TR
        heartTrigPoints
        respTrigPoints
        respPercentage = 30
        respWindow
        
        % Rates
        heartRateTime
        heartRateTimeFiltered
        respRateTime
        respRateTimeFiltered
        meanHeartRate
        meanRespRate
        
    end
    
    
    % ---------------------------------------------------------------------------------
    % Public methods
    % -------------------------------------------------------------------------------
    
    methods (Access = public)
        
        
        
        % ---------------------------------------------------------------------------------
        % Object constructor
        % ---------------------------------------------------------------------------------
        function obj = retroNav
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Read navigator up/down flip switch value
        % ---------------------------------------------------------------------------------
        function obj = readFlipSwitch(obj,flipswitchvalue)
            
            switch flipswitchvalue
                
                case 'Up'
                    obj.up_down = 1;
                    
                case 'Down'
                    obj.up_down = -1;
                    
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Extract the navigator signals
        % ---------------------------------------------------------------------------------
        function objNav = extractNavigator(objNav,objData)
            
            switch objData.dataType
                
                case '2D'
                    extractNavigator2D;
                    
                case {'3D','3Dp2roud'}
                    extractNavigator3D;
                    
                case '2Dms'
                    extractNavigator2Dms;
                    
                case '2Dradial'
                    extractNavigatorRadial;
                    
            end
            
            
            % ---------------------------------------------------------------------------------
            % ----- 2D single-slice data -------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function amplitude = extractNavigator2D
                
                objNav.nav_amplitude = cell(objData.nr_coils);
                
                for coilnr = 1:objData.nr_coils
                    
                    % extracts the navigator data from the raw k-space data
                    % outputs a 1D array of doubles
                    
                    % size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nr_eperiments,dimz,dimy,dimx] = size(objData.data{coilnr});
                    
                    % extract the navigator and put it in a long array
                    navdata_amplitude = reshape(permute(objData.data{coilnr},[3,2,1,4]),nr_eperiments*dimy*dimz,dimx);
                    
                    if objData.nr_nav_points_used > 1
                        
                        % take the principal component of the data
                        data = navdata_amplitude(:,objData.primary_navigator_point-objData.nr_nav_points_used+1:objData.primary_navigator_point);
                        [coeff,~,~] = pca(data);
                        dataPCA = data*coeff;
                        
                        % take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';
                        
                    else
                        
                        % single nav point
                        amplitude = abs(navdata_amplitude(:,objData.primary_navigator_point))';
                        
                    end
                    
                    % detrend
                    amplitude(1,:) = detrend(amplitude(1,:));
                    
                    % make a guess whether the respiration peaks are positive or negative
                    nr_elements = length(amplitude);
                    first_element = round(0.4*nr_elements);
                    last_element = round(0.8*nr_elements);
                    max_amplitude = abs(max(detrend(amplitude(1,first_element:last_element))));
                    min_amplitude = abs(min(detrend(amplitude(1,first_element:last_element))));
                    if min_amplitude > max_amplitude
                        amplitude(1,:) = -amplitude(1,:);
                    end
                    
                    % multiple with +1 or -1 depending on switch
                    amplitude = amplitude * objNav.up_down;
                    
                    % return the final nav amplitude array
                    objNav.nav_amplitude{coilnr} = amplitude;
                    
                end
                
            end % extractNavigator2D
            
            
            
            
            % ---------------------------------------------------------------------------------
            % ----- 2D multi-slice data -------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function amplitude = extractNavigator2Dms
                
                objNav.nav_amplitude = cell(objData.nr_coils);
                
                for coilcounter = 1:objData.nr_coils
                    
                    % extracts the navigator data from the raw k-space data of multi-slice 2D data
                    % outputs a 1D array of doubles
                    
                    % size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nr_experiments,dimz,dimy,dimx] = size(objData.data{coilcounter});
                    
                    % extract the navigator and put it in a long array
                    % y-dimension, repetitions, slice, readout
                    navdata_amplitude = reshape(permute(objData.data{coilcounter},[3,1,2,4]),nr_experiments*dimy*dimz,dimx);
                    
                    if objData.nr_nav_points_used > 1
                        
                        % Take the principal component of the data
                        data = navdata_amplitude(:,objData.primary_navigator_point-objData.nr_nav_points_used+1:objData.primary_navigator_point);
                        [coeff,~,~] = pca(data);
                        dataPCA = data*coeff;
                        
                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';
                        
                    else
                        
                        % single nav point
                        amplitude = abs(navdata_amplitude(:,objData.primary_navigator_point))';
                        
                    end
                    
                    % detrend
                    amplitude(1,:) = detrend(amplitude(1,:));
                    
                    % make a guess whether the respiration peaks are positive or negative in the different slices
                    % this will not be needed with out-of-slice navigator
                    
                    nr_elements = length(amplitude);
                    nr_elementsperslice = round(nr_elements/dimz);
                    
                    for i = 1:dimz
                        
                        % First and last navigator point for each slice
                        first_element = (i-1)*nr_elementsperslice + 1;
                        last_element = i*nr_elementsperslice;
                        
                        % Only look at part of that data away from the start to prevent transient effects
                        first_element1 = (i-1)*nr_elementsperslice + 1 + round(0.4*nr_elementsperslice);
                        last_element1 = i*nr_elementsperslice - round(0.1*nr_elementsperslice);
                        
                        % Min/max of navigator
                        max_amplitude = abs(max(detrend(amplitude(1,first_element1:last_element1))));
                        min_amplitude = abs(min(detrend(amplitude(1,first_element1:last_element1))));
                        
                        if min_amplitude > max_amplitude
                            amplitude(1,first_element:last_element) = -amplitude(1,first_element:last_element);
                        end
                        
                        amplitude(1,first_element:last_element) = detrend(amplitude(1,first_element:last_element));
                        
                    end
                    
                    % multiple with +1 or -1 depending on global switch
                    amplitude = amplitude * objNav.up_down;
                    
                    % return the final nav amplitude
                    objNav.nav_amplitude{coilcounter} = amplitude;
                    
                end
                
            end % extractNavigator2Dms
            
            
            % ---------------------------------------------------------------------------------
            % ----- 3D data -------------------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function amplitude = extractNavigator3D
                
                objNav.nav_amplitude = cell(objData.nr_coils);
                
                for coilnr = 1:objData.nr_coils
                    
                    % extracts the navigator data from the raw k-space data
                    % outputs a 1D array of doubles
                    
                    % size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nr_rep,dimz,dimy,dimx] = size(objData.data{coilnr});
                    
                    % extract the navigator and put it in a long array
                    navdata_amplitude = reshape(permute(objData.data{coilnr},[2,3,1,4]),nr_rep*dimy*dimz,dimx);
                    
                    if objData.nr_nav_points_used > 1
                        
                        % take the principal component of the data
                        data = navdata_amplitude(:,objData.primary_navigator_point-objData.nr_nav_points_used+1:objData.primary_navigator_point);
                        [coeff,~,~] = pca(data);
                        dataPCA = data*coeff;
                        
                        % take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';
                        
                    else
                        
                        % single nav point
                        amplitude = abs(navdata_amplitude(:,objData.primary_navigator_point))';
                        
                    end
                    
                    % detrend
                    amplitude(1,:) = detrend(amplitude(1,:));
                    
                    % make a guess whether the respiration peaks are positive or negative
                    nr_elements = length(amplitude);
                    first_element = round(0.4*nr_elements);
                    last_element = round(0.6*nr_elements);
                    max_amplitude = abs(max(amplitude(1,first_element:last_element)));
                    min_amplitude = abs(min(amplitude(1,first_element:last_element)));
                    if min_amplitude > max_amplitude
                        amplitude = -amplitude;
                    end
                    
                    % multiple with +1 or -1 depending on switch
                    amplitude = amplitude * objNav.up_down;
                    
                    % return the final nav amplitude
                    objNav.nav_amplitude{coilnr} = amplitude;
                    
                end
                
            end % extractNavigator3D
            
            
            
            % ---------------------------------------------------------------------------------
            % ----- Radial data ---------------------------------------------------------------
            % ---------------------------------------------------------------------------------
            function extractNavigatorRadial
                
                objNav.nav_amplitude = cell(objData.nr_coils);
                
                for coilnr = 1:objData.nr_coils
                    
                    % extracts the navigator data from the raw k-space data
                    % outputs a 1D array of doubles
                    
                    % size of the input data, 4th dimension is the readout direction which contains the navigator
                    [nr_rep,dimz,dimy,dimx] = size(objData.data{coilnr});
                    
                    % determine the phase offset of the individual spokes based on the center navigator point
                    phaseoffset = zeros(nr_rep*dimz*dimy,2);
                    cnt1 = 1;
                    for i = 1:nr_rep
                        for j = 1:dimz
                            for k = 1:dimy
                                phaseoffset(cnt1,:) = [k,angle(objData.data{coilnr}(i,j,k,objData.primary_navigator_point))];
                                cnt1 = cnt1 + 1;
                            end
                        end
                    end
                    
                    %phaseoffset = sort(phaseoffset,1);
                    %phaseoffset = [unique(phaseoffset(:,1)),accumarray(phaseoffset(:,1),phaseoffset(:,2),[],@mean)];
                    
                    %phaseoffset_x = phaseoffset(:,1);
                    %phaseoffset_y = phaseoffset(:,2);
                    %phaseoffset_y = unwrap(phaseoffset_y);
                    %figure(1)
                    %scatter(phaseoffset_x,phaseoffset_y)
                    
                    % extract the navigator and put it in a long array
                    navdata_amplitude = reshape(permute(objData.data{coilnr},[3,2,1,4]),nr_rep*dimy*dimz,dimx);
                    
                    if objData.nr_nav_points_used > 1
                        
                        range = round(objData.nr_nav_points_used/2);
                        
                        % Take the principal component of the data
                        data = navdata_amplitude(:,objData.primary_navigator_point-range:objData.primary_navigator_point+range);
                        [coeff,~,~] = pca(data);
                        dataPCA = data*coeff;
                        
                        % Take the principal component of the data
                        amplitude = abs(dataPCA(:,1))';
                        phase = angle(dataPCA(:,1))';
                        
                    else
                        
                        % single nav point
                        amplitude = abs(navdata_amplitude(:,objData.primary_navigator_point))';
                        phase = angle(navdata_amplitude(:,objData.primary_navigator_point))';
                        
                    end
                    
                    % Detrend
                    amplitude(1,:) = detrend(amplitude(1,:));
                    
                    % Make a guess whether the respiration peaks are positive or negative
                    nr_elements = length(amplitude);
                    first_element = round(0.4*nr_elements);
                    last_element = round(0.6*nr_elements);
                    max_amplitude = abs(max(amplitude(1,first_element:last_element)));
                    min_amplitude = abs(min(amplitude(1,first_element:last_element)));
                    if min_amplitude > max_amplitude
                        amplitude = -amplitude;
                    end
                    
                    % multiple with +1 or -1 depending on switch
                    amplitude = amplitude * objNav.up_down;
                    
                    % return the final nav amplitude
                    objNav.nav_amplitude{coilnr} = amplitude;
                    objNav.nav_phase{coilnr} = phase;
                    
                end
                
            end % extractNavigatorRadial
            
            
        end % extractNavigator
       
       
       
        % ---------------------------------------------------------------------------------
        % Determine the power-frequency spectrum of the navigator
        % ---------------------------------------------------------------------------------
        function  objNav = determinePowerSpectrumPCA(objNav,objData)
            
            if objData.nr_coils > 1
                
                data = zeros([length(objNav.nav_amplitude{1}),objData.nr_coils]);
                for i = 1:objData.nr_coils
                    data(:,i) = objNav.nav_amplitude{i};
                end
                
                % take the principal component of the data
                [coeff,~,~] = pca(data);
                dataPCA = data*coeff;
                amplitude = dataPCA(:,1);
                
            else
                
                amplitude = objNav.nav_amplitude{1}';
                
            end
            
            % include only those navigators in when includewindow == 1 and excludewindow == 1
            amplitude = amplitude.*objData.includeWindow.*objData.excludeWindow;
            
            % determine the frequency power spectrum
            y = fft(amplitude);
            fs = 1000/objData.TR;                       % sample frequency in Hz
            n = length(amplitude);                      % number of samples
            objNav.frequency = (0:n-1)*(fs/n)*60;       % frequency range in bpm
            power = abs(y).^2/n;                        % power of the DFT
            
            % determine frequency and harmonics of k-space trajectory and set those to zero
            kfreq = 60/(0.001*objData.NO_VIEWS*objData.TR);
            ifreq = (fs/n)*60;
            
            for i = 1:10
                power(round(i*kfreq/ifreq)+1) = 0;
                power(round(i*kfreq/ifreq))   = 0;
                power(round(i*kfreq/ifreq)-1) = 0;
            end
            
            % smooth the power spectrum with moving average
            power = movmean(power,3);
            objNav.powerspectrum = power;
            
            % detect heart rate
            minheartbpm = objNav.physiofiltersettings(1);
            maxheartbpm = objNav.physiofiltersettings(2);
            minidx = round(minheartbpm*n/(fs*60));
            maxidx = round(maxheartbpm*n/(fs*60));
            [~, idx] = max(power(minidx:maxidx));
            objNav.detectedHR = round(idx*fs*60/n + minheartbpm);
            
            % detect respiratory rate
            minRRbpm = objNav.physiofiltersettings(3);
            maxRRbpm = objNav.physiofiltersettings(4);
            minidx = round(minRRbpm*n/(fs*60));
            maxidx = round(maxRRbpm*n/(fs*60));
            [~, idx] = max(power(minidx:maxidx));
            objNav.detectedRR = round(idx*fs*60/n + minRRbpm);
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Filter the navigator
        % ---------------------------------------------------------------------------------
        function objNav = filterNavPCA(objNav,objData)
            
            % applies a bandwidth filter on the navigator data
            
            sf = 1000/objData.TR;                   % sampling frequency in Hz = 1/TR[ms]
            resp_harmonics = 2;                     % number of higher order harmonics for respiratory frequency, 2 = main + 1 harmonic
            order = objNav.physiofiltersettings(5); % butterworth filter order
            
            if objData.nr_coils > 1
                
                data = zeros([length(objNav.nav_amplitude{1}),objData.nr_coils]);
                for i = 1:objData.nr_coils
                    data(:,i) = objNav.nav_amplitude{i};
                end
                
                % take the principal component of the data
                [coeff,~,~] = pca(data);
                dataPCA = data*coeff;
                
                % take the principal component of the data
                amplitude = dataPCA(:,1);
                
            else
                
                amplitude = objNav.nav_amplitude{1}';
                
            end
            
            % filter for heart motion
            hrf = objNav.detectedHR/60;       % expected heartrate in Hz = hr[bpm]/60
            bwh = objNav.bandwidthHR/60;      % bandwidth heartrate in Hz = [bpm]/60
            [b,a] = butter(order,[hrf-0.5*bwh,hrf+0.5*bwh]/(sf/2),'bandpass');      % butterworth bandpass filter
            heart_outputdata = filtfilt(b,a,amplitude);                             % apply zero-phase shift filtering
            
            % detrend
            heart_outputdata = detrend(heart_outputdata);
            
            % normalize envelope
            factor = round(100/hrf);      % adapt the envelope setting to match the expected heart rate frequency
            [env,~] = envelope(heart_outputdata,factor,'peak');
            objNav.heartNav = heart_outputdata./abs(env);
            
            % filter for respiration motion
            while true
                
                rrf=objNav.detectedRR/60;       % expected resprate in Hz = rr[bpm]/60
                bwr=objNav.bandwidthRR/60;      % bandwidth resprate in Hz = [bpm]/60
                
                resp_outputdata = zeros(size(amplitude));
                
                if objNav.detectedRR<45
                    [b, a] = butter(order,(rrf+0.5*bwr)/(sf/2),'low');           % butterworth lowpass filter for low frequencies
                    resp_outputdata = filtfilt(b,a,amplitude);
                else
                    for i = 1:resp_harmonics
                        [b, a] = butter(order,[i*rrf-0.5*bwr,i*rrf+0.5*bwr]/(sf/2),'bandpass');       % butterworth bandpass filter
                        resp_outputdata = resp_outputdata + (1/i^2)*filtfilt(b,a,amplitude);          % apply zero-phase shift filtering
                    end
                end
                
                % In some cases the filter produces NaN when filtering with too low
                % frequency. In those cases the respirate and bandwidth will be
                % increased until there are no more NaN
                
                if sum(isnan(resp_outputdata))==0
                    break;
                end
                
                objNav.detectedRR = objNav.detectedRR+1;
                objNav.bandwidthRR = objNav.bandwidthRR+1;
                
            end
            
            % detrend and normalize envelope
            resp_outputdata = detrend(resp_outputdata);
            
            % normalize envelope
            factor = round(150/rrf);  % adapt the envelope setting to match the expected respiration rate frequency
            [env,~] = envelope(resp_outputdata,factor,'peak');
            objNav.respNav = resp_outputdata./abs(env);
            
            % return the principal component, or in case 1 coil the orignal navigator data, and detrend
            objNav.PCANav = detrend(amplitude);
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Heart filter shape for plotting
        % ---------------------------------------------------------------------------------
        function output = filterShape(objNav,objData,order,type)
            
            % determines the filter shape for display purposes
            % order = order of the filter
            % type = 'heart' or 'resp'
            
            sf = 1000/objData.TR;           % sampling frequency in Hz = 1/TR[ms]
            
            if strcmp(type,'resp')
                hr = objNav.detectedRR;
                bw = objNav.bandwidthRR;
            else
                hr = objNav.detectedHR;
                bw = objNav.bandwidthHR;
            end
            
            cf=hr/60;      % expected heartrate in Hz = hr[bpm]/60
            bw1=bw/60;      % bandwidth in Hz = [bpm]/60
            
            % check that bw1 is smaller than center-frequency
            if cf<=bw1
                bw1 = 0.98*cf;
            end
            
            if hr < 45
                [b, a] = butter(order,(cf+0.5*bw1)/(sf/2),'low');           % butterworth lowpass filter for low frequencies
            else
                [b,a] = butter(order,[cf-0.5*bw1,cf+0.5*bw1]/(sf/2),'bandpass');   % butterworth bandpass filter
            end
            
            [h,w] = freqz(b,a,2000,sf);
            output = [w*60, abs(h)];
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Determine the cardiac trigger points
        % ---------------------------------------------------------------------------------
        function objNav = trigPointsHeart(objNav,objData)
            
            % extracts the ECG trigger points from the navigators
            
            % find the peaks and locations = fast
            objNav.heartTrigPoints = retroNav.peakfinder(objNav.heartNav',[],[],[],false,true);
            
            % trigger points are in units of samples (actual time is thus heartTrigPoints*TR)
         
            % backup plan in case peakfinder fails = slower
            if length(objNav.heartTrigPoints)<20
                
                % minimal distance 50% of expected heart rate [in points]
                if(isnan(objNav.meanHeartRate))
                    objNav.meanHeartRate = 500;
                end
                dist = 0.50*(60/objNav.meanHeartRate)/(objData.TR/1000);
                interpolationfactor = 60;
                
                % cubic spline interpolation of the data
                nrl = length(objNav.heartNav);
                navi = interp1(1:nrl,objNav.heartNav(1:nrl),1:1/interpolationfactor:nrl,'spline');
                
                % find the peaks and locations
                [~,locs]=findpeaks(navi,'MinPeakDistance',dist*interpolationfactor);
                locs = locs + interpolationfactor/2;
                
                % recalculate orginal fractional time point peak positions
                objNav.heartTrigPoints = locs/interpolationfactor;
                
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Determine the respiration trigger points
        % ---------------------------------------------------------------------------------
        function objNav = trigPointsResp(objNav)
            
            % extracts the ECG trigger points from the navigators
            
            % find the peaks and locations
            objNav.respTrigPoints = retroNav.peakfinder(objNav.respNav',[],[],[],false,true);
            
            % trigger points are in units of samples (actual time is thus heartTrigPoints*TR)
            
            % backup plan in case peakfinder fails
            if length(objNav.respTrigPoints)<10
                
                interpolationfactor = objNav.splineFactor;
                nrl = size(objNav.respNav,1);
                navi = interp1(1:nrl,objNav.respNav,1:1/interpolationfactor:nrl,'spline');
                
                % find the peaks and locations
                [~,locs] = findpeaks(navi,'MinPeakProminence',0.1);
                locs = locs + interpolationfactor/2;
                
                % recalculate orginal fractional time point peak positions
                objNav.respTrigPoints = locs/interpolationfactor;
                
            end
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Determine average heart rate
        % ---------------------------------------------------------------------------------
        function objNav = calcHeartRate(objNav,objData)
            
            % determine heart rate as function of time in bpm
            
            nrlocs = size(objNav.heartTrigPoints,2);
            
            rate = zeros(nrlocs-1,1);
            for i=1:nrlocs-1
                rate(i) = 60/((objNav.heartTrigPoints(i+1)-objNav.heartTrigPoints(i))*(objData.TR/1000));
            end
            
            hr = rate';
            hrf = movmedian(hr,32); % smooth trendline
            
            % determine the mean cardiac rate
            includedata = objData.includeWindow.*objData.excludeWindow;                     % data-window which is included
            includedata = round(resample(includedata,4*length(hr),length(includedata)));    % resample the data-window to nr samples heartrate
            includedata = round(resample(includedata,length(hr),length(includedata)));      % in 2 steps, to prevent an overflow error
            
            objNav.meanHeartRate = round(median(nonzeros(includedata.*hr)));                % take the median of the heartrate
            objNav.heartRateTime = hr;
            objNav.heartRateTimeFiltered = hrf;
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Determine average respiration rate
        % ---------------------------------------------------------------------------------
        function objNav = calcRespRate(objNav,objData)
            
            % determine respiration rate as function of time in bpm
            
            nrlocs = size(objNav.respTrigPoints,2);
            
            rate = zeros(nrlocs-1,1);
            for i=1:nrlocs-1
                rate(i) = 60/((objNav.respTrigPoints(i+1)-objNav.respTrigPoints(i))*(objData.TR/1000));
            end
            
            resp = rate';
            respf = movmedian(resp,32);     % smooth trendline
            
            % determine the mean respiration rate
            includedata = objData.includeWindow.*objData.excludeWindow;                     % data-window which is included
            includedata = round(resample(includedata,4*length(resp),length(includedata)));  % resample the data-window to nr samples respiration rate
            includedata = round(resample(includedata,length(resp),length(includedata)));    % in 2 steps, to prevent an overflow error
            
            objNav.meanRespRate = round(median(nonzeros(includedata.*resp)));               % take the median of the respirationrate
            objNav.respRateTime = resp;
            objNav.respRateTimeFiltered = respf;
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Determine respiration windows
        % ---------------------------------------------------------------------------------
        function objNav = makeRespWindow(objNav,objData)
            
            %(rlocs,mean_resp,resp_percentage,tr,nr_klines)
            
            rlocs = objNav.respTrigPoints;
            mean_resp = mean(objNav.respRateTimeFiltered);
            resp_percentage = objNav.respPercentage;
            tr = objData.TR;
            nr_klines = objData.nrKlines;
            
            % app.retroNavPars.respTrigPoints,mean(app.retroNavPars.respRateTimeFiltered),app.app.retroNavPars.respPercentage,app.retroPars.TR,app.nr_klines
            
            % this function creates an array (time line) of rectangular boxes of 0's and 1's around the detected respiratory signals
            % later on 1 means that there is a respiration, for which the k-lines will be discarded
            
            resp_a = round(0.5*(resp_percentage/100)*(60/mean_resp)*1000/tr);   % 1/2 window width around respiratory peak locations
            
            window = zeros(nr_klines,1);   % fill array with zeros
            
            for i = 1 : size(rlocs,2)
                center = round(rlocs(1,i));
                for j = center - resp_a : center + resp_a
                    if (j>0) && (j<=nr_klines)
                        window(j)=1;
                    end
                end
            end
            
            objNav.respWindow = window;
            
        end
        
        
        
    end
        
    
    
    
    % ---------------------------------------------------------------------------------
    % Static methods
    % ---------------------------------------------------------------------------------
    
    methods (Static)
        
        function varargout = peakfinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)
            %PEAKFINDER Noise tolerant fast peak finding algorithm
            %   INPUTS:
            %       x0 - A real vector from the maxima will be found (required)
            %       sel - The amount above surrounding data for a peak to be,
            %           identified (default = (max(x0)-min(x0))/4). Larger values mean
            %           the algorithm is more selective in finding peaks.
            %       thresh - A threshold value which peaks must be larger than to be
            %           maxima or smaller than to be minima.
            %       extrema - 1 if maxima are desired, -1 if minima are desired
            %           (default = maxima, 1)
            %       includeEndpoints - If true the endpoints will be included as
            %           possible extrema otherwise they will not be included
            %           (default = true)
            %       interpolate - If true quadratic interpolation will be performed
            %           around each extrema to estimate the magnitude and the
            %           position of the peak in terms of fractional indicies. Note that
            %           unlike the rest of this function interpolation assumes the
            %           input is equally spaced. To recover the x_values of the input
            %           rather than the fractional indicies you can do:
            %           peakX = x0 + (peakLoc - 1) * dx
            %           where x0 is the first x value and dx is the spacing of the
            %           vector. Output peakMag to recover interpolated magnitudes.
            %           See example 2 for more information.
            %           (default = false)
            %
            %   OUTPUTS:
            %       peakLoc - The indicies of the identified peaks in x0
            %       peakMag - The magnitude of the identified peaks
            
            % Perform error checking and set defaults if not passed in
            narginchk(1, 6);
            nargoutchk(0, 2);
            
            s = size(x0);
            flipData =  s(1) < s(2);
            len0 = numel(x0);
            if len0 ~= s(1) && len0 ~= s(2)
                error('PEAKFINDER:Input','The input data must be a vector')
            elseif isempty(x0)
                varargout = {[],[]};
                return;
            end
            if ~isreal(x0)
                warning('PEAKFINDER:NotReal','Absolute value of data will be used')
                x0 = abs(x0);
            end
            
            if nargin < 2 || isempty(sel)
                sel = (max(x0)-min(x0))/4;
            elseif ~isnumeric(sel) || ~isreal(sel)
                sel = (max(x0)-min(x0))/4;
                warning('PEAKFINDER:InvalidSel',...
                    'The selectivity must be a real scalar.  A selectivity of %.4g will be used',sel)
            elseif numel(sel) > 1
                warning('PEAKFINDER:InvalidSel',...
                    'The selectivity must be a scalar.  The first selectivity value in the vector will be used.')
                sel = sel(1);
            end
            
            if nargin < 3 || isempty(thresh)
                thresh = [];
            elseif ~isnumeric(thresh) || ~isreal(thresh)
                thresh = [];
                warning('PEAKFINDER:InvalidThreshold',...
                    'The threshold must be a real scalar. No threshold will be used.')
            elseif numel(thresh) > 1
                thresh = thresh(1);
                warning('PEAKFINDER:InvalidThreshold',...
                    'The threshold must be a scalar.  The first threshold value in the vector will be used.')
            end
            
            if nargin < 4 || isempty(extrema)
                extrema = 1;
            else
                extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
                if extrema == 0
                    error('PEAKFINDER:ZeroMaxima','Either 1 (for maxima) or -1 (for minima) must be input for extrema');
                end
            end
            
            if nargin < 5 || isempty(includeEndpoints)
                includeEndpoints = true;
            end
            
            if nargin < 6 || isempty(interpolate)
                interpolate = false;
            end
            
            x0 = extrema*x0(:); % Make it so we are finding maxima regardless
            thresh = thresh*extrema; % Adjust threshold according to extrema.
            dx0 = diff(x0); % Find derivative
            dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
            ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign
            
            % Include endpoints in potential peaks and valleys as desired
            if includeEndpoints
                x = [x0(1);x0(ind);x0(end)];
                ind = [1;ind;len0];
                minMag = min(x);
                leftMin = minMag;
            else
                x = x0(ind);
                minMag = min(x);
                leftMin = min(x(1), x0(1));
            end
            
            % x only has the peaks, valleys, and possibly endpoints
            len = numel(x);
            
            if len > 2 % Function with peaks and valleys
                % Set initial parameters for loop
                tempMag = minMag;
                foundPeak = false;
                
                if includeEndpoints
                    signDx = sign(diff(x(1:3)));
                    if signDx(1) <= 0 % The first point is larger or equal to the second
                        if signDx(1) == signDx(2) % Want alternating signs
                            x(2) = [];
                            ind(2) = [];
                            len = len-1;
                        end
                    else % First point is smaller than the second
                        if signDx(1) == signDx(2) % Want alternating signs
                            x(1) = [];
                            ind(1) = [];
                            len = len-1;
                        end
                    end
                end
                
                % Skip the first point if it is smaller so we always start on a maxima
                if x(1) >= x(2)
                    ii = 0;
                else
                    ii = 1;
                end
                
                % Preallocate max number of maxima
                maxPeaks = ceil(len/2);
                peakLoc = zeros(maxPeaks,1);
                peakMag = zeros(maxPeaks,1);
                cInd = 1;
                % Loop through extrema which should be peaks and then valleys
                while ii < len
                    ii = ii+1; % This is a peak
                    % Reset peak finding if we had a peak and the next peak is bigger
                    %   than the last or the left min was small enough to reset.
                    if foundPeak
                        tempMag = minMag;
                        foundPeak = false;
                    end
                    
                    % Found new peak that was lager than temp mag and selectivity larger
                    %   than the minimum to its left.
                    if x(ii) > tempMag && x(ii) > leftMin + sel
                        tempLoc = ii;
                        tempMag = x(ii);
                    end
                    
                    % Make sure we don't iterate past the length of our vector
                    if ii == len
                        break; % We assign the last point differently out of the loop
                    end
                    
                    ii = ii+1; % Move onto the valley
                    % Come down at least sel from peak
                    if ~foundPeak && tempMag > sel + x(ii)
                        foundPeak = true; % We have found a peak
                        leftMin = x(ii);
                        peakLoc(cInd) = tempLoc; % Add peak to index
                        peakMag(cInd) = tempMag;
                        cInd = cInd+1;
                    elseif x(ii) < leftMin % New left minima
                        leftMin = x(ii);
                    end
                end
                
                % Check end point
                if includeEndpoints
                    if x(end) > tempMag && x(end) > leftMin + sel
                        peakLoc(cInd) = len;
                        peakMag(cInd) = x(end);
                        cInd = cInd + 1;
                    elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
                        peakLoc(cInd) = tempLoc;
                        peakMag(cInd) = tempMag;
                        cInd = cInd + 1;
                    end
                elseif ~foundPeak
                    if x(end) > tempMag && x(end) > leftMin + sel
                        peakLoc(cInd) = len;
                        peakMag(cInd) = x(end);
                        cInd = cInd + 1;
                    elseif tempMag > min(x0(end), x(end)) + sel
                        peakLoc(cInd) = tempLoc;
                        peakMag(cInd) = tempMag;
                        cInd = cInd + 1;
                    end
                end
                
                % Create output
                if cInd > 1
                    peakInds = ind(peakLoc(1:cInd-1));
                    peakMags = peakMag(1:cInd-1);
                else
                    peakInds = [];
                    peakMags = [];
                end
            else % This is a monotone function where an endpoint is the only peak
                [peakMags,xInd] = max(x);
                if includeEndpoints && peakMags > minMag + sel
                    peakInds = ind(xInd);
                else
                    peakMags = [];
                    peakInds = [];
                end
            end
            
            % Apply threshold value.  Since always finding maxima it will always be
            %   larger than the thresh.
            if ~isempty(thresh)
                m = peakMags>thresh;
                peakInds = peakInds(m);
                peakMags = peakMags(m);
            end
            
            if interpolate && ~isempty(peakMags)
                middleMask = (peakInds > 1) & (peakInds < len0);
                noEnds = peakInds(middleMask);
                
                magDiff = x0(noEnds + 1) - x0(noEnds - 1);
                magSum = x0(noEnds - 1) + x0(noEnds + 1)  - 2 * x0(noEnds);
                magRatio = magDiff ./ magSum;
                
                peakInds(middleMask) = peakInds(middleMask) - magRatio/2;
                peakMags(middleMask) = peakMags(middleMask) - magRatio .* magDiff/8;
            end
            
            % Rotate data if needed
            if flipData
                peakMags = peakMags.';
                peakInds = peakInds.';
            end
            
            % Change sign of data if was finding minima
            if extrema < 0
                peakMags = -peakMags;
                x0 = -x0;
            end
            
            % Plot if no output desired
            if nargout == 0
                if isempty(peakInds)
                    disp('No significant peaks found')
                else
                    figure;
                    plot(1:len0,x0,'.-',peakInds,peakMags,'ro','linewidth',2);
                end
            else
                varargout = {peakInds,peakMags};
            end
            
        end % peakfinder
        
        
    end
    
    
end