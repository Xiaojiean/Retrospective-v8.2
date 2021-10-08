function newrpr = makerpr(inputrpr,par)


parameters = {
    ':EDITTEXT LAST_ECHO ',':EDITTEXT MAX_ECHO ', ...
    ':EDITTEXT LAST_EXPT ',':EDITTEXT MAX_EXPT ', ...
    ':EDITTEXT SAMPLES_DIM1 ',':EDITTEXT DATA_LENGTH1 ', ':EDITTEXT OUTPUT_SIZE1 ', ...
    ':EDITTEXT SAMPLES_DIM2 ',':EDITTEXT DATA_LENGTH2 ', ':EDITTEXT OUTPUT_SIZE2 ', ...
    ':EDITTEXT SAMPLES_DIM3 ',':EDITTEXT DATA_LENGTH3 ', ':EDITTEXT OUTPUT_SIZE3 ', ...
    ':EDITTEXT LAST_SLICE ',':EDITTEXT MAX_SLICE ', ...
    ':COMBOBOX FFT_DIM1 ',':COMBOBOX FFT_DIM2 ',':COMBOBOX FFT_DIM3 ', ...
    ':RADIOBUTTON VIEW_ORDER_2'
    };
            
replacepars = {par.NoEchoes,par.NoEchoes, ...
    par.NoExperiments, par.NoExperiments, ...
    par.NoSamples, par.NoSamples, par.NoSamples, ...
    par.NoViews, par.NoViews, par.NoViews, ...
    par.NoViews2, par.NoViews2, par.NoViews2, ...
    par.NoSlices, par.NoSlices, ...
    par.NoSamples, par.NoViews, par.NoViews2, ...
    par.View2order
    };
            
            


for i = 1:length(parameters)
    
    txt = parameters{i};
    var = replacepars{i};
            
    pos = strfind(inputrpr,txt);
    
    if ~isempty(pos)
        
        if ~isstring(var)
            
            oldtxtlength = strfind(inputrpr(pos+length(txt):pos+length(txt)+15),char(13))-1;
            newtext = [num2str(var),'     '];
            newtext = newtext(1:6);
            inputrpr = replaceBetween(inputrpr,pos+length(txt),pos+length(txt)+oldtxtlength-1,newtext);
            
        else
            
            oldtxtlength = strfind(inputrpr(pos+length(txt):pos+length(txt)+15),char(13))-1;
            newtext = strcat(" ",var,"           ");
            newtext = extractBefore(newtext,12);
            inputrpr = replaceBetween(inputrpr,pos+length(txt),pos+length(txt)+oldtxtlength-1,newtext);
            
        end
         
    end
    
end

newrpr = inputrpr;

end