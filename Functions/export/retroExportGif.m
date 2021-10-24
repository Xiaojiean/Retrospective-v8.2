function folder_name = retroExportGif(parameters,gifExportPath,movie,tag,window,level,recoType,acqDur)
% Exports movie to animated gif


% Dimensions
[nrFrames,~,~,nrSlices,nrDynamics] = size(movie);


% Create folder if not exist, and clear
folder_name = [gifExportPath,[filesep,'RETRO_GIFS_',num2str(nrFrames),'_',num2str(nrSlices),'_',num2str(nrDynamics),'_',tag]];
if (~exist(folder_name, 'dir')); mkdir(folder_name); end
delete([folder_name,filesep,'*']);


% Scale from 0 to 255
window = window*255/max(movie(:));
level = level*255/max(movie(:));
movie = movie*255/max(movie(:));


% Window and level
movie = (255/window)*(movie - level + window/2);
movie(movie < 0) = 0;
movie(movie > 255) = 255;


% Correct for non-square aspect ratio
gifImageSize = 512; % size of longest axis
dimy = gifImageSize;
dimx = round(dimy * parameters.aspectratio);
if parameters.PHASE_ORIENTATION
    dimx = gifImageSize;
    dimy = round(dimx * parameters.aspectratio);
end
fct = max([dimx dimy]);
dimx = round(gifImageSize * dimx / fct);
dimy = round(gifImageSize * dimy / fct);


% Variable flip-angle
if parameters.VFA_size > 1
    dynamiclabel = '_flipangle_';
else
    dynamiclabel = '_dynamic_';
end


if strcmp(recoType,'realtime')

    % Dynamic movie

    delayTime = acqDur/nrDynamics;

    for i = 1:nrSlices

        slice = ['0',num2str(i)];
        slice = slice(end-1:end);

        for j = 1:nrFrames

            dyn = ['00',num2str(j)];
            dyn = dyn(end-2:end);

            for idx = 1:nrDynamics

                image = uint8(squeeze(movie(j,:,:,i,idx)));
                image = imresize(image,[dimx,dimy]);

                if idx == 1
                    imwrite(image,[folder_name,filesep,'movie_',tag,'_slice_',slice,'frame',dyn,'.gif'],'DelayTime',delayTime,'LoopCount',inf);
                else
                    imwrite(image,[folder_name,filesep,'movie_',tag,'_slice_',slice,'frame',dyn,'.gif'],'WriteMode','append','DelayTime',delayTime);
                end

            end

        end

    end

else

    % Frames movie

    delayTime = 1/nrFrames;

    for i = 1:nrSlices

        slice = ['0',num2str(i)];
        slice = slice(end-1:end);

        for j = 1:nrDynamics

            dyn = ['00',num2str(j)];
            dyn = dyn(end-2:end);

            for idx = 1:nrFrames

                imaget = uint8(squeeze(movie(idx,:,:,i,j)));
                image(i,:,:,idx,j) = imresize(imaget,[dimx,dimy]); %#ok<AGROW> 
                imagegif = squeeze(image(i,:,:,idx,j));

                if idx == 1
                    imwrite(imagegif,[folder_name,filesep,'movie_',tag,'_slice_',slice,dynamiclabel,dyn,'.gif'],'DelayTime',delayTime,'LoopCount',inf);
                else
                    imwrite(imagegif,[folder_name,filesep,'movie_',tag,'_slice_',slice,dynamiclabel,dyn,'.gif'],'WriteMode','append','DelayTime',delayTime);
                end

            end

        end

    end

    if nrSlices > 1

        % Make image collage
        dv = divisor(nrSlices);
        rows = dv(round(end/2));
        cols = nrSlices/rows;
        imageCollection = zeros(rows*dimx,cols*dimy,nrFrames,nrDynamics);
        for i = 1:rows
            for j = 1:cols
                imageCollection((i-1)*dimx+1:i*dimx,(j-1)*dimy+1:j*dimy,:,:) = squeeze(image((j-1)*i+i,:,:,:,:));
            end
        end

        % Export collage
        for j = 1:nrDynamics

            dyn = ['00',num2str(j)];
            dyn = dyn(end-2:end);

            for idx = 1:nrFrames

                imagegif = squeeze(imageCollection(:,:,idx,j));

                if idx == 1
                    imwrite(imagegif,[folder_name,filesep,'collage_',tag,dynamiclabel,dyn,'.gif'],'DelayTime',delayTime,'LoopCount',inf);
                else
                    imwrite(imagegif,[folder_name,filesep,'collage_',tag,dynamiclabel,dyn,'.gif'],'WriteMode','append','DelayTime',delayTime);
                end

            end
        end

    end

end



    function d = divisor(n)
        % Provides a list of integer divisors of a number

        % Find prime factors of number
        pf = factor(n);         % Prime factors of n
        upf = unique(pf);       % Unique

        % Calculate the divisors
        d = upf(1).^(0:1:sum(pf == upf(1)))';
        for f = upf(2:end)
            d = d*(f.^(0:1:sum(pf == f)));
            d = d(:);
        end
        d = sort(d)';  

    end



end