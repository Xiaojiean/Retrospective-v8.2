classdef retroReco
    
    % Reconstruction and movie class for retrospective app
    
    properties
        
        movieExp                        % movie for movie export
        movieApp                        % movie for viewing in the app
        senseMap                        % sense map
        sos = true                      % Sum of squares reco true or false
        multiSliceFlag = false          % multi-slice true or false
        multiDynamicFlag = false        % mutli-dynamic true or false
        
    end
    
    
    
    
    
    methods
        
        % ---------------------------------------------------------------------------------
        % Object constructor
        % ---------------------------------------------------------------------------------
        function obj = retroReco
            
        end % retroReco
        
        
        
        % ---------------------------------------------------------------------------------
        % 2D reconstruction
        % ---------------------------------------------------------------------------------
        function objReco = reco2D(objReco, objData, objKspace, app)
            
            if app.bartDetectedFlag
                
                % Perform CS reconstruction with the Bart toolbox, preferred option
                app.TextMessage('Reconstructing the data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;
                cs_reco2D_mc;
                app.ProgressGauge.Value = 100;
                drawnow;
                
            else
                
                % Perform CS reconstruction in matlab, slow but works
                app.TextMessage('WARNING: Reconstructing the data with MATLAB ...');
                app.TextMessage('Slow reconstruction, please be patient ...');
                app.ProgressGauge.Value = 0;
                app.SetStatus(1);
                drawnow;
                cs_reco2D_mat_mc;
                app.ProgressGauge.Value = 100;
                drawnow;
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 2D compressed sensing reconstruction with Matlab
            % ---------------------------------------------------------------------------------
            function cs_reco2D_mat_mc
                
                kspace_in = objKspace.kSpace;
                averages = objKspace.kSpaceAvg;
                lambda_TV = app.TVcineEditField.Value;
                
                % kspace = {coil}[CINE, y, x, slice, dynamic]
                %                   1   2  3    4       5
                nr_cine = size(kspace_in{1},1);
                nr_slices = size(kspace_in{1},4);
                nr_dynamics = size(kspace_in{1},5);
                nc = objData.nr_coils;
                
                % in case of 1 frame, duplicate that frame to facilitate reconstruction
                if nr_cine == 1
                    for i = 1:nc
                        kspace_in{i}(2,:,:,:,:) = kspace_in{i}(1,:,:,:,:);
                    end
                end
                
                % kspace data y,x,t,slices,dynamics
                for i = 1:nc
                    kspace_in{i} = permute(kspace_in{i},[2,3,1,4,5]);
                end
                
                % kspace data y,x,t,slices,dynamics,coils
                kspace = zeros([size(kspace_in{1}),nc]);
                for i = 1:nc
                    kspace(:,:,:,:,:,i) = kspace_in{i};
                end
                
                % averages data y,x,t,slices,dynamics
                averages = permute(averages,[2,3,1,4,5]);
                
                % reset progress counter
                param.iteration = 0;
                
                % pad to next power of 2
                dimy = 2^nextpow2(size(kspace,1));
                dimx = 2^nextpow2(size(kspace,2));
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                dimf = size(kspace,3);
                
                % pre-allocate memory for image_out
                image_out = zeros(dimf,dimy,dimx,nr_slices,nr_dynamics);
                
                % slice and dynamic loop
                for slice = 1:nr_slices
                    
                    for dynamic = 1:nr_dynamics
                        
                        % kspace of slice and dynamic
                        kdata = squeeze(kspace(:,:,:,slice,dynamic,:));
                        
                        % zero padding
                        padsizey = round((dimy - size(kdata,1))/2);
                        padsizex = round((dimx - size(kdata,2))/2);
                        kdata = padarray(kdata,[padsizey,padsizex,0],'both');
                        kdata = kdata(1:dimy,1:dimx,:,:);
                        
                        % size of the data
                        [ny,nx,~,nc]=size(kdata);
                        
                        % normalize the data in the range of approx 0 - 1 for better numerical stability
                        kdata = kdata/max(abs(kdata(:)));
                        
                        % kspace mask, 0 = nodata, 1 = data, zero-pad to same size as k-space
                        mask = squeeze(averages(:,:,:,slice,dynamic));
                        mask = padarray(mask,[padsizey,padsizex,0],'both');
                        mask = mask(1:dimy,1:dimx,:);
                        mask = mask./mask;
                        mask(isnan(mask)) = 1;
                        mask = logical(mask);
                        
                        % coil sensitivity map
                        b1 = ones(ny,nx,nc);
                        
                        % data
                        param.y = kdata;
                        
                        % reconstruction design matrix
                        param.E = Emat_yxt(mask,b1);
                        
                        % Total variation (TV) constraint in the temporal domain
                        % for 'consistency' with Bart reconstruction, TV seems to be scale empirically by a factor of 8
                        % TV only in the time domain
                        param.TV = TVOP;
                        param.TVWeight = lambda_TV/8;
                        
                        % number of iterations, 2 x 10 iterations
                        param.nite = 10;
                        param.nouter = 2;
                        param.totaliterations = nr_slices * nr_dynamics * param.nouter * param.nite;
                        
                        % linear reconstruction
                        kdata1 = randn(size(kdata))/2000 + kdata;  % add a little bit of randomness, such that linear reco is not exactly right
                        recon_dft = param.E'*kdata1;
                        
                        % iterative reconstruction
                        recon_cs = recon_dft;
                        for n = 1:param.nouter
                            [recon_cs,param.iteration] = CSL1NlCg(app,recon_cs,param);
                        end
                        
                        % rearrange to correct orientation: frames, x, y
                        image_tmp = flip(permute(squeeze(abs(recon_cs)),[3, 1, 2]),2);
                        
                        % output reconstructed image
                        image_out(:,:,:,slice,dynamic) = image_tmp;
                        
                    end
                    
                end
                
                % correct back to 1 frame reconstruction
                if nr_cine == 1
                    image_out = image_out(1,:,:,:,:);
                end
                
                % shift image in phase-encoding direction if needed
                objReco.movieExp = circshift(image_out,objData.pixelshift1,2);
                objReco.senseMap = ones(size(objReco.movieExp));
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 2D compressed sensing reconstruction with Bart toolbox
            % ---------------------------------------------------------------------------------
            function cs_reco2D_mc
                
                kspace_in = objKspace.kSpace;
                Wavelet = app.WVxyzEditField.Value;
                TVxy = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                ESPIRiT = app.ESPIRiTCheckBox.Value;
                SOS = objReco.sos;
                
                % app = matlab app
                % kspace_in = sorted k-space
                % nc = number of RF receiver coils
                % averages_in = number of averages per k-line
                % Wavelet = wavelet L1-norm regularization factor
                % TVt = total variation in CINE dimension
                % TVxy = total variation in xy-dimension regularization
                % TVd = total variation in dynamic dimension
                % ESPIRiT = reconstruction of multi-coil data with ESPIRiT true or false
                
                % kspace_in = {coil}[CINE, x, y, slice, dynamic]
                %                     1    2  3    4       5
                dimc = size(kspace_in{1},1);
                dimx = 2^nextpow2(size(kspace_in{1},2));
                dimy = 2^nextpow2(size(kspace_in{1},3));
                dimz = size(kspace_in{1},4);
                dimd = size(kspace_in{1},5);
                nc = objData.nr_coils;
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                
                % resize k-space to next power of 2
                for i = 1:nc
                    kspace_in{i} = bart(app,['resize -c 1 ',num2str(dimx),' 2 ',num2str(dimy)],kspace_in{i});
                end
                
                % put all data in a normal matrix
                kspace = zeros(dimc,dimx,dimy,dimz,dimd);
                for l = 1:dimd
                    for k = 1:dimz
                        for j = 1:dimc
                            for i = 1:nc
                                kspace(j,:,:,k,l,i) = kspace_in{i}(j,:,:,k,l);
                            end
                        end
                    end
                end
                
                % Bart dimensions  Bart   Matlab
                % 	READ_DIM,       0       1   z
                % 	PHS1_DIM,       1       2   y
                % 	PHS2_DIM,       2       3   x
                % 	COIL_DIM,       3       4   coils
                % 	MAPS_DIM,       4       5   sense maps
                % 	TE_DIM,         5       6
                % 	COEFF_DIM,      6       7
                % 	COEFF2_DIM,     7       8
                % 	ITER_DIM,       8       9
                % 	CSHIFT_DIM,     9       10
                % 	TIME_DIM,       10      11  cardiac / respiratory CINE frames
                % 	TIME2_DIM,      11      12  dynamics
                % 	LEVEL_DIM,      12      13
                % 	SLICE_DIM,      13      14  slices
                % 	AVG_DIM,        14      15
                
                kspace_pics = permute(kspace,[7,3,2,6,8,9,10,11,12,13,1,5,14,4]);
                
                % wavelet in spatial dimensions 2^1+2^2=6
                % total variation in spatial dimensions 2^1+2^2=6
                % total variation in cine dimension 2^10 = 1024
                % total variation in dynamic dimension 2^11 = 2048
                
                if ESPIRiT && nc>1
                    
                    % ESPIRiT reconstruction
                    TextMessage(app,'ESPIRiT reconstruction ...');
                    
                    % Calculate coil sensitivity maps with ecalib bart function
                    kspace_pics_sum = sum(kspace_pics,[11,12]);
                    sensitivities = bart(app,'ecalib -S -I -a', kspace_pics_sum);      % ecalib with softsense
                    
                    picscommand = 'pics ';
                    if Wavelet>0
                        picscommand = [picscommand,' -RW:6:0:',num2str(Wavelet)];
                    end
                    if TVxy>0
                        picscommand = [picscommand,' -RT:6:0:',num2str(TVxy)];
                    end
                    if LR>0
                        % Locally low-rank in the spatial domain
                        blocksize = round(dimx/16);  % Block size
                        picscommand = [picscommand,' -RL:6:6:',num2str(LR),' -b',num2str(blocksize)];
                    end
                    if TVt>0
                        picscommand = [picscommand,' -RT:1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picscommand = [picscommand,' -RT:2048:0:',num2str(TVd)];
                    end
                    image_reg = bart(app,picscommand,kspace_pics,sensitivities);
                    
                    % Sum of squares reconstruction
                    if SOS
                        image_reg = bart(app,'rss 16', image_reg);
                    end
                    image_reg = abs(image_reg);
                    
                else
                    
                    % Reconstruction without sensitivity correction, including coil scaling
                    sensitivities = ones(1,dimy,dimx,nc,1,1,1,1,1,1,1,1,1,dimz);
                    
                    % regular reconstruction
                    picscommand = 'pics ';
                    if Wavelet>0
                        picscommand = [picscommand,' -RW:6:0:',num2str(Wavelet)];
                    end
                    if TVxy>0
                        picscommand = [picscommand,' -RT:6:0:',num2str(TVxy)];
                    end
                    if LR>0
                        picscommand = [picscommand,' -b20 -RL:6:6:',num2str(LR)];
                    end
                    if TVt>0
                        picscommand = [picscommand,' -RT:1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picscommand = [picscommand,' -RT:2048:0:',num2str(TVd)];
                    end
                    image_reg = bart(app,picscommand,kspace_pics,sensitivities);
                    
                    % Take absolute values
                    image_reg = abs(image_reg);
                    
                end
                
                % rearrange to correct orientation: frames, x, y, slices, dynamics, coils
                image_reg = reshape(image_reg,[dimy,dimx,dimc,dimd,dimz]);
                image_out = flip(permute(image_reg,[3,2,1,5,4]),3);
                
                % sense map orientations: x, y, slice, map1, map2
                sensemap = flip(permute(abs(sensitivities),[3,2,14,4,5,1,6,7,8,9,10,11,12,13]),2);
                
                % normalize sense map to reasonable value range
                sensemap = sensemap*4095/max(sensemap(:));
                
                % shift image in phase-encoding direction if needed
                objReco.movieExp = circshift(image_out,objData.pixelshift1,2);
                objReco.senseMap = circshift(sensemap,objData.pixelshift1,2);
                
            end
            
        end
        
                
        
        % ---------------------------------------------------------------------------------
        % 3D reconstruction
        % ---------------------------------------------------------------------------------
        function objReco = reco3D(objReco, objData, objKspace, app)
            
            if app.bartDetectedFlag
                
                % perform CS reconstruction with the Bart toolbox, preferred option
                app.TextMessage('Reconstructing the 3D data with the BART toolbox ...');
                app.ProgressGauge.Value = 0;
                drawnow;
                cs_reco3D_mc;
                app.ProgressGauge.Value = 100;
                drawnow;
                
            else
                
                % perform cs reconstruction in matlab, slow but works
                app.TextMessage('WARNING: Reconstructing the data with MATLAB ...');
                app.TextMessage('Slow reconstruction, please be patient ...');
                app.ProgressGauge.Value = 0;
                app.SetStatus(1);
                cs_reco3D_mat_mc;
                app.ProgressGauge.Value = 100;
                drawnow;
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 3D compressed sensing reconstruction with Matlab
            % ---------------------------------------------------------------------------------
            function cs_reco3D_mat_mc
                
                kspace_in = objKspace.kSpace;
                averages = objKspace.kSpaceAvg;
                lambda_TV = app.TVcineEditField.Value;
                
                % kdata_in = {coil}[CINE, y, x, z, dynamic]
                %                     1   2  3  4     5
                nr_cine = size(kspace_in{1},1);
                nr_dynamics = size(kspace_in{1},5);
                nc = objData.nr_coils;
                
                % in case of 1 frame, duplicate that frame to facilitate reconstruction
                if nr_cine == 1
                    for i = 1:nc
                        kspace_in{i}(2,:,:,:,:) = kspace_in{i}(1,:,:,:,:);
                    end
                end
                
                % kspace data z,y,x,cine,dynamic
                for i = 1:nc
                    kspace_in{i} = permute(kspace_in{i},[4,2,3,1,5]);
                end
                
                % kspace data z,y,x,cine,dynamic,coils
                kspace = zeros([size(kspace_in{1}),nc]);
                for i = 1:nc
                    kspace(:,:,:,:,:,i) = kspace_in{i};
                end
                
                % averages data z,y,x,cine,dynamics
                averages = permute(averages,[4,2,3,1,5]);
                
                % reset progress counter
                param.iteration = 0;
                
                % pad to next power of 2
                dimz = 2^nextpow2(size(kspace,1));
                dimy = 2^nextpow2(size(kspace,2));
                dimx = 2^nextpow2(size(kspace,3));
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                dimf = size(kspace,4);
                
                % pre-allocate memory for image_out
                image_out = zeros(dimf,dimx,dimy,dimz,nr_dynamics);
                
                for dynamic = 1:nr_dynamics
                    
                    % kspace of slice and dynamic
                    kdata = squeeze(kspace(:,:,:,:,dynamic,:));
                    
                    % zero padding
                    padsizez = round((dimz - size(kdata,1))/2);
                    padsizey = round((dimy - size(kdata,2))/2);
                    padsizex = round((dimx - size(kdata,3))/2);
                    kdata = padarray(kdata,[padsizez,padsizey,padsizex,0],'both');
                    kdata = kdata(1:dimz,1:dimy,1:dimx,:,:);
                    
                    % size of the data, z,y,x,frames,coils
                    [nz,ny,nx,~,nc] = size(kdata);
                    
                    % normalize the data in the range of approx 0 - 1 for better numerical stability
                    kdata = kdata/max(abs(kdata(:)));
                    
                    % kspace mask, 0 = nodata, 1 = data, zero-pad to same size as k-space
                    mask = squeeze(averages(:,:,:,:,dynamic));
                    mask = padarray(mask,[padsizez,padsizey,padsizex,0],'both');
                    mask = mask(1:dimz,1:dimy,1:dimx,:);
                    mask = mask./mask;
                    mask(isnan(mask)) = 1;
                    mask = logical(mask);
                    
                    % coil sensitivity map
                    b1 = ones(nz,ny,nx,nc);
                    
                    % data
                    param.y = kdata;
                    
                    % reconstruction design matrix
                    param.E = Emat_zyxt(mask,b1);
                    
                    % Total variation (TV) constraint in the temporal domain
                    % for 'consistency' with Bart reconstruction, TV seems to be scale empirically by a factor of 8
                    % TV only in the time domain
                    param.TV = TVOP3D;
                    param.TVWeight = lambda_TV/8;
                    
                    % number of iterations, 2 x 10 iterations
                    param.nite = 10;
                    param.nouter = 2;
                    param.totaliterations = nr_dynamics * param.nouter * param.nite;
                    
                    % linear reconstruction
                    kdata1 = randn(size(kdata))/2000 + kdata;  % add a little bit of randomness, such that linear reco is not exactly right
                    recon_dft = param.E'*kdata1;
                    
                    % iterative reconstruction
                    recon_cs=recon_dft;
                    for n=1:param.nouter
                        [recon_cs,param.iteration] = CSL1NlCg(app,recon_cs,param);
                    end
                    
                    % rearrange to correct orientation: frames, x, y, z
                    image_tmp = flip(permute(abs(recon_cs),[4, 2, 3, 1]),2);
                    
                    % output reconstructed image
                    image_out(:,:,:,:,dynamic) = image_tmp;
                    
                end
                
                % correct back to 1 frame reconstruction
                if nr_cine == 1
                    image_out = image_out(1,:,:,:,:);
                end
                
                % shift image in phase-encoding direction if needed
                objReco.movieExp = circshift(image_out,objData.pixelshift1,2);
                objReco.movieExp = circshift(image_out,objData.pixelshift2,1);
                objReco.senseMap = ones(size(objReco.movieExp));
                
            end
            
            
            % ---------------------------------------------------------------------------------
            % 3D compressed sensing reconstruction with Bart toolbox
            % ---------------------------------------------------------------------------------
            function cs_reco3D_mc
                
                kspace_in = objKspace.kSpace;
                Wavelet = app.WVxyzEditField.Value;
                TVxyz = app.TVxyzEditField.Value;
                LR = app.LLRxyzEditField.Value;
                TVt = app.TVcineEditField.Value;
                TVd = app.TVdynEditField.Value;
                ESPIRiT = app.ESPIRiTCheckBox.Value;
                
                % app = matlab app
                % kspace_in = sorted k-space
                % nc = number of RF receiver coils
                % averages_in = number of averages per k-line
                % Wavelet = wavelet L1-norm regularization factor
                % TVt = total variation in time regularization
                % TVxyz = total variation in xyz-dimension regularization
                % VNSA = variable number of signal averages correction true or false
                % ESPIRiT = reconstruction of multi-coil data with ESPIRiT true or false
                
                % kspace_in = {coil}[CINE, x, y, slice, dynamic]
                %                     1    2  3    4       5
                dimc = size(kspace_in{1},1);
                dimx = 2^nextpow2(size(kspace_in{1},2));
                dimy = 2^nextpow2(size(kspace_in{1},3));
                dimz = 2^nextpow2(size(kspace_in{1},4));
                dimd = size(kspace_in{1},5);
                nc = objData.nr_coils;
                
                % for convenience make rectangular matrix size
                mdimxy = max([dimx,dimy]);
                dimx = mdimxy;
                dimy = mdimxy;
                
                % resize to next power of 2
                for i = 1:nc
                    kspace_in{i} = bart(app,['resize -c 1 ',num2str(dimx),' 2 ',num2str(dimy),' 3 ',num2str(dimz)],kspace_in{i});
                end
                
                % kspace suitable for bart
                kspace = zeros(dimc,dimx,dimy,dimz,dimd,nc);
                for l = 1:dimd
                    for k = 1:dimz
                        for j = 1:dimc
                            for i = 1:nc
                                kspace(j,:,:,k,l,i) = kspace_in{i}(j,:,:,k,l);
                            end
                        end
                    end
                end
                
                % Bart dimensions
                % 	READ_DIM,       1   z
                % 	PHS1_DIM,       2   y
                % 	PHS2_DIM,       3   x
                % 	COIL_DIM,       4   coils
                % 	MAPS_DIM,       5   sense maps
                % 	TE_DIM,         6
                % 	COEFF_DIM,      7
                % 	COEFF2_DIM,     8
                % 	ITER_DIM,       9
                % 	CSHIFT_DIM,     10
                % 	TIME_DIM,       11  cardiac / respiratory CINE
                % 	TIME2_DIM,      12  dynamics
                % 	LEVEL_DIM,      13
                % 	SLICE_DIM,      14  slices
                % 	AVG_DIM,        15
                
                kspace_pics = permute(kspace,[4,3,2,6,7,8,9,10,11,12,1,5,13,14]);
                
                % wavelet in spatial dimensions 2^0+2^1+2^2=7
                % total variation in spatial dimensions 2^0+2^1+2^2=7
                % total variation in time 2^10 = 1024
                % total variation in dynamic dimension 2^11 = 2048
                
                if ESPIRiT && nc>1
                    
                    TextMessage(app,'ESPIRiT reconstruction ...');
                    kspace_pics_sum = sum(kspace_pics,[11,12]);
                    sensitivities = bart(app,'ecalib -I -S -a', kspace_pics_sum);
                    
                    app.ProgressGauge.Value = 25;
                    drawnow;
                    
                    picscommand = 'pics -S ';
                    if Wavelet>0
                        picscommand = [picscommand,' -RW:7:0:',num2str(Wavelet)];
                    end
                    if TVxyz>0
                        picscommand = [picscommand,' -RT:7:0:',num2str(TVxyz)];
                    end
                    if LR>0
                        % Locally low-rank in the spatial domain
                        blocksize = round(dimx/16);  % Block size
                        picscommand = [picscommand,' -RL:7:7:',num2str(LR),' -b',num2str(blocksize)];
                    end
                    if TVt>0
                        picscommand = [picscommand,' -RT:1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picscommand = [picscommand,' -RT:2048:0:',num2str(TVd)];
                    end
                    image_reg = bart(app,picscommand,kspace_pics,sensitivities);
                    
                    app.ProgressGauge.Value = 95;
                    drawnow;
                    image_reg = bart(app,'rss 16', image_reg);
                    image_reg = abs(image_reg);
                    
                else
                    
                    % Reconstruction without sensitivity correction
                    sensitivities = ones(dimz,dimy,dimx,nc,1,1,1,1,1,1,1,1,1,1);
                    
                    picscommand = 'pics -S ';
                    if Wavelet>0
                        picscommand = [picscommand,' -RW:7:0:',num2str(Wavelet)];
                    end
                    if TVxyz>0
                        picscommand = [picscommand,' -RT:7:0:',num2str(TVxyz)];
                    end
                    if LR>0
                        picscommand = [picscommand,' -RL:7:7:',num2str(LR)];
                    end
                    if TVt>0
                        picscommand = [picscommand,' -RT:1024:0:',num2str(TVt)];
                    end
                    if TVd>0
                        picscommand = [picscommand,' -RT:2048:0:',num2str(TVd)];
                    end
                    
                    image_reg = abs(bart(app,picscommand,kspace_pics,sensitivities));
                    
                end
                
                % rearrange to correct orientation: frames, x, y, z, dynamics
                image_reg = reshape(image_reg,[dimz,dimy,dimx,dimc,dimd]);
                image_out = flip(permute(image_reg,[4,3,2,1,5]),3);
                
                % sense map orientations: x, y, z, map1, map2
                sensemap = flip(permute(abs(sensitivities),[3,2,1,4,5,6,7,8,9,10,11,12,13,14]),2);
                
                % normalize sense map to reasonable value range
                sensemap = sensemap*4095/max(sensemap(:));
                
                % shift image in phase-encoding direction if needed
                objReco.movieExp = circshift(image_out,objData.pixelshift1,2);
                objReco.movieExp = circshift(image_out,objData.pixelshift2,1);
                objReco.senseMap = circshift(sensemap,objData.pixelshift1,2);
                
            end
            
        end % reco3D
        
                
        
        % ---------------------------------------------------------------------------------
        % Radial reconstruction with the Bart toolbox
        % ---------------------------------------------------------------------------------
        function objReco = recoRadial(objReco, objData, objKspace, app)
            
            app.TextMessage('Reconstructing 2D radial data with the BART toolbox ...');
            app.ProgressGauge.Value = 0;
            drawnow;
            
            kspace_in = objKspace.kSpace;
            Wavelet = app.WVxyzEditField.Value;
            TVxy = app.TVxyzEditField.Value;
            TVt = app.TVcineEditField.Value;
            nc = objData.nr_coils;
            dimx = objData.dimx;
            dimy = objData.dimy;
            ksteps = objData.nrKsteps; %#ok<NASGU> 
            
            % app = matlab app
            % kspace_in = sorted k-space
            % nc = number of RF receiver coils
            % averages_in = number of averages per k-line
            % Wavelet = wavelet L1-norm regularization factor
            % TVt = total variation in time regularization
            % TVxy = total variation in xy-dimension regularization
            
            % radial k-space trajectory for Bart toolbox
            radial = 360;           % 360 degrees
            u_spokes = true;        % unique spokes true/false
            half_spoke = false;     %#ok<NASGU> % half trajectory true/false
            
            for n = 1:dimy
                
                % caluculate angle
                if radial == 180
                    a(n) = (n-1)*180/dimy; %#ok<*AGROW> 
                elseif ~u_spokes
                    a(n) = (n-1)*360/dimy;
                elseif (n-1) < dimy/2
                    a(n) = (n-1)*360/dimy;
                else
                    a(n) = (n-1)*360/dimy + 180/dimy;
                end
                
                % calculate x and y values
                for m = 1:dimx
                    r(m,n) = (dimx-1)/2 -(m-1);
                    x = r(m,n)*-cos(a(n)*(pi/180));
                    y = r(m,n)*sin(a(n)*(pi/180));
                    traj(1,m,n) = x;
                    traj(2,m,n) = y;
                    traj(3,m,n) = 0;
                end
                
            end
            
            % bart (kz,ky,kx,channels,multiple sense, ..., ..., cardiac frames at pos 11)
            % sense maps for 2D data:  (kz,ky,kx,....)
            for i = 1:nc
                kspace(:,:,:,:,i) = kspace_in{i};
            end
            
            %figure(1)
            %imshow(real(squeeze(kspace(1,:,:))),[-2*pi,2*pi]);
            %imshow(abs(squeeze(kspace(1,:,:))),[0 10]);
            
            % rearrange for the bart toolbox
            kspace_pics = permute(kspace,[4,3,2,5,6,7,8,9,10,11,1]);
            
            % create trajectory
            %bartcommand = ['traj -r -y',num2str(ksteps),' -x',num2str(dimx),' -q0:0:0'];
            %traj = bart(bartcommand);
            
            % sensitivity map
            kspace_pics_sum = sum(kspace_pics,11);
            lowres_img = bart(app,'nufft -i -l6 -d16:16:1 -t', traj, kspace_pics_sum);
            lowres_ksp = bart(app,'fft -u 7', lowres_img);
            
            highres_img = abs(bart(app,'nufft -i -t -l6 -d128:128:1 -t', traj, kspace_pics_sum)); %#ok<NASGU> 
            figure(1)
            imshow(abs(lowres_img),[0,1.5*max(abs(lowres_img(:)))]);
            
            % zeropad to full size
            bartcommand = ['resize -c 0 ',num2str(dimx),' 1 ',num2str(dimx)];
            ksp_zerop = bart(app, bartcommand, lowres_ksp);
            
            % calculate sensitivity map
            sensitivities = bart(app, 'ecalib -m1', ksp_zerop); %#ok<NASGU> 
            
            %disp(size(sensitivities))
            %sensitivities = ones(128,128);
            %disp(size(sensitivities))
            
            sensitivities = ones(dimx,dimx,1,1,1,1,1,1,1,1,1,1,1);
            
            picscommand = ['pics -S -u1 -RW:6:0:',num2str(Wavelet),' -RT:6:0:',num2str(TVxy),' -RT:1024:0:',num2str(TVt),' -t'];
            image_reg = bart(app,picscommand,traj,kspace_pics,sensitivities); %#ok<NASGU> 
            
            image_reg = bart(app,'nufft -i -l1 -d128:128:1 -t', traj, kspace_pics);
            
            % rearrange to correct orientation: frames,x, y
            image_reg = abs(squeeze(image_reg(:,:,1,1,1,1,1,1,1,1,:)));
            image_out = flip(permute(image_reg,[3, 2, 1]),3);
            sensemap = flip(permute(squeeze(abs(sensitivities)),[2, 1, 3, 4]),2);
            
            % normalize sense map to reasonable value range
            sensemap = sensemap*1024/max(sensemap(:));
            
            % shift image in phase-encoding direction if needed
            objReco.movieExp = circshift(image_out,objData.pixelshift1,2);
            objReco.senseMap = circshift(sensemap,objData.pixelshift1,2);
            
            app.ProgressGauge.Value = 100;
            drawnow;
            
        end
        
        
        
        % ---------------------------------------------------------------------------------
        % Normalize movie intensity
        % ---------------------------------------------------------------------------------
        function objReco = normImages(objReco)
            
            % normalize the images to 2^15 range
            
            objReco.movieExp = round(32766*objReco.movieExp/max(objReco.movieExp(:)));
            objReco.movieApp = objReco.movieExp;
            
        end % normImages
        
        

        % ---------------------------------------------------------------------------------
        % Normalize movie intensity
        % ---------------------------------------------------------------------------------
        function objReco = determineMultiDimensions(objReco)
            
            % multislice
            if size(objReco.movieApp,4) > 1
                objReco.multiSliceFlag = true;
            else
                objReco.multiSliceFlag = false;
            end

            % multidynamic
            if size(objReco.movieApp,5) > 1
                objReco.multiDynamicFlag = true;
            else
                objReco.multiDynamicFlag = false;
            end

        end % determineMultiDimensions
        

        
        % ---------------------------------------------------------------------------------
        % Left-right or up-down phase orientation
        % ---------------------------------------------------------------------------------
        function objReco = phaseOrientation(objReco, objData, app)
            
            % Phase_orientation correction for viewing
            
            if objData.PHASE_ORIENTATION
                app.TextMessage('INFO: phase orientation = 1 / horizontal ... ')
                objReco.movieApp = permute(rot90(permute(objReco.movieApp,[2,3,4,1,5]),1),[4,1,2,3,5]);
                objReco.movieApp = flip(objReco.movieApp,3);
                objReco.senseMap = rot90(objReco.senseMap,1);
                objReco.senseMap = flip(objReco.senseMap,2);
            else
                app.TextMessage('INFO: phase orientation = 0 / vertical ... ')
            end
            
        end % phaseOrientation
        
        

        % ---------------------------------------------------------------------------------
        % Image reconstruction: SUR files
        % ---------------------------------------------------------------------------------
        function obj = recoSurFiles(obj, surpath, suffix, mrdfilename, rprfilename)

            % SUR file names
            surfiles = [surpath, suffix, '_00###.SUR'];

            % Link with the server
            m_Recon = actxserver('recon.Application');

            set(m_Recon,'Visible',1);
            set(m_Recon,'DisplayImages',1);

            % Filenames
            set(m_Recon,'DataFile',mrdfilename);
            set(m_Recon,'RPRFile',rprfilename);
            set(m_Recon,'ImageFile',surfiles);

            % Delete old SUR files
            scmd = ['del /Q ', surpath, '*.SUR'];
            system(scmd);

            % Do the reco
            invoke(m_Recon,'Run');

            % Wait for recon to complete
            while get(m_Recon,'StatusID')~=4
                if get(m_Recon,'StatusID')==5
                    break;
                end
                pause(0.1);
            end

            % Stop the link
            invoke(m_Recon,'Quit');

        end % recoSurFiles







    end % methods
    
end % retroReco

