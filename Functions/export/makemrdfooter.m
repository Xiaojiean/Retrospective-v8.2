function footer = makemrdfooter(inputFooter,par)


parameters = {':NO_SAMPLES no_samples, ',':NO_VIEWS no_views, ',':NO_VIEWS_2 no_views_2, ', ... 
                ':NO_ECHOES no_echoes, ',':EXPERIMENT_ARRAY no_experiments, ',':NO_AVERAGES no_averages, ', ...
                ':VAR pe1_order, ',':VAR slice_nav, ',':VAR radial_on, ', ... 
                ':VAR frame_loop_on, ',':VAR tr, ',':VAR te, ', ...
                ':BATCH_SLICES batch_slices, ',':NO_SLICES no_slices, ', ...
                ':VAR VFA_size, ',':VAR ti, ',':VAR pe2_centric_on, ', ...
                ':VAR tr_extra_us, '
                };

replacePars = {par.NoSamples,par.NoViews,par.NoViews2, ... 
                par.NoEchoes,par.NoExperiments,par.NoAverages, ... 
                par.peorder,par.slicenav,par.radialon, ... 
                par.frameloopon,par.tr,par.te, ...
                par.batchslices,par.NoSlices, ...
                par.vfasize, par.ti, par.pe2_centric_on, ...
                par.tr_extra_us
                };

            
% Replace all simple valued parameters
for i = 1:length(parameters)
    
    txt = parameters{i};
    var = replacePars{i};
    
    pos = strfind(inputFooter,txt);
    
    if ~isempty(pos)
        oldTextLength = strfind(inputFooter(pos+length(txt):pos+length(txt)+6),char(13))-1;
        newText = [num2str(var),'     '];
        newText = newText(1:6);
        inputFooter = replaceBetween(inputFooter,pos+length(txt),pos+length(txt)+oldTextLength-1,newText);
    end
    
end

% VFA array replace
newText = '';
for i = 1:length(par.vfaangles)
    newText = newText + ", " + num2str(par.vfaangles(i));
    if mod(i,8) == 2
       newText = newText + newline; 
    end
end
txt = ':VAR_ARRAY VFA_angles, ';
pos1 = strfind(inputFooter,txt);
newStr = extractAfter(inputFooter,pos1);
pos2 = strfind(newStr,':');
pos3 = pos1+pos2(1)-2;
inputFooter = replaceBetween(inputFooter,pos1+length(txt)-2,pos3,newText);


% Replace the footer with the new footer
footer = inputFooter;


end