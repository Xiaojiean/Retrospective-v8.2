
function struct = jcampread(filename)

% Based on: 
%   readreco.m
%   Function to read  RECO files to a structure


%   Check if filename and showimg where given as argument, show UI if filename was not.
if nargin == 0
    [fname,pname] = uigetfile({'*.*', 'All Files (*.*)'},'Select JCAMP DX file', 100,100);
    filename = [pname fname];
end

% Open file read-only big-endian
[fid,message]=fopen(filename,'r','b');
skipline=0;

% Loop through separate lines
if fid~=-1
    while 1
        if skipline
            line=nextline;
            skipline=0;
        else
            line=fgetl(fid);
        end
        % Testing the text lines 
        while length(line)<2
            line=fgetl(fid);
        end
        %% Parameters and optional size of parameter are on lines starting with '##'
        if line(1:2) == '##'
            %% Parameter extracting and formatting
            % Read parameter name
            paramname = fliplr(strtok(fliplr(strtok(line,'=')),'#'));
            % Check for illegal parameter names starting with '$' and correct (Matlab does not accepts variable names starting with $)
            if paramname(1) == '$'
                paramname = paramname(2:length(paramname));
                % Check if EOF, if true return    
            elseif paramname(1:3) == 'END'
                break
            end
            %% Parameter value formatting
            paramvalue = fliplr(strtok(fliplr(line),'='));
           
            %paramname
            % Check if parameter values are in a matrix and read the next line
            if paramvalue(1) == '('
                paramvaluesize = str2num(fliplr(strtok(fliplr(strtok(paramvalue,')')),'(')));
                % Create an empty matrix with size 'paramvaluesize' check if only one dimension
                if ~isempty(paramvaluesize)
                    if size(paramvaluesize,2) == 1
                        paramvaluesize = [paramvaluesize,1];
                    end
                    % Read the next line
                    nextline = fgetl(fid);
                    % See whether next line contains a character array
                    if nextline(1) == '<' 
                        paramvalue = fliplr(strtok(fliplr(strtok(nextline,'>')),'<'));
                    elseif strcmp(nextline(1),'L') || strcmp(nextline(1),'A') || strcmp(nextline(1),'H')
                        paramvalue = nextline;
                    else
                        % Check if matrix has more then one dimension
                        if paramvaluesize(2) ~= 1
                            paramvaluelong = str2num(nextline);
                            while (length(paramvaluelong)<(paramvaluesize(1)*paramvaluesize(2))) & (nextline(1:2) ~= '##')
                                nextline = fgetl(fid);
                                paramvaluelong = [paramvaluelong str2num(nextline)];
                            end
                            if (length(paramvaluelong)==(paramvaluesize(1)*paramvaluesize(2))) & (~isempty(paramvaluelong))
                                paramvalue=reshape(paramvaluelong,paramvaluesize(1),paramvaluesize(2));
                            else
                                paramvalue=paramvaluelong;
                            end
                            if length(nextline)>1
                                if (nextline(1:2) ~= '##') 
                                    skipline=1; 
                                end
                            end
                        else
                            % If only 1 dimension just assign whole line to paramvalue
                            paramvalue = str2num(nextline);
                            if ~isempty(str2num(nextline))
                                while length(paramvalue)<paramvaluesize(1)
                                    line=fgetl(fid);
                                    paramvalue = [paramvalue str2num(line)];
                                end
                            end
                        end
                    end
                else
                    paramvalue='';
                end
            end
            
            % Add paramvalue to structure.paramname
            if isempty(findstr(paramname,'_'))
                eval(['struct.' paramname '= paramvalue;']);
              else  
                try
                    eval(['struct.' lower(paramname(1:findstr(paramname,'_')-1)) '.' lower(paramname(findstr(paramname,'_')+1:length(paramname))) '= paramvalue;']);
                catch
                    eval(['struct.' lower(paramname(1:findstr(paramname,'_')-1)) '.' datestr(str2num(paramname(findstr(paramname,'_')+1:findstr(paramname,'_')+2)),9) ...
                            paramname(findstr(paramname,'_')+2:length(paramname)) '= paramvalue;']);
                end
            end   
        elseif line(1:2) == '$$'
            % The two $$ lines are not parsed for now
        end
    end
    % Close file
    fclose(fid);
else
    disp(message)
    struct=-1;
end
