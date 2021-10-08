% Description: Function to convert multidimensional complex data to MRD format file
% Author: Ruslan Garipov / MR Solutions Ltd
% Date: 17/04/2020
% Inputs: string filename, N-dimensional data matrix, dimensions structure with the following fields:
% .NoExperiments
% .NoEchoes
% .NoSlices
% .NoSamples
% .NoViews
% .NoViews2
% footer - int8 data type footer as copied from an MRD containing a copy of
% the PPR file, including preceeding 120-byte zeros
% Output: 1 if write was successful, 0 if not, -1 if failed early (e.g. the dimension checks)

function writedata2mrd(filename, data, dimensions, footer)
% data: multidimensional, complex float, with dimensions arranged
% dimensions: structure

% get dimensions of the actual image data 
    if (size(data,1)~=dimensions.NoSamples)
        return;
    end
    if (size(data,2)>1)
        if (size(data,2)~=dimensions.NoViews)
            return;
        end
    end
    if (size(data,3)>1)
        if (size(data,3)~=dimensions.NoViews2)
            return;
        end
    end
    if (size(data,4)>1)
        if (size(data,4)~=dimensions.NoSlices)
            return;
        end
    end
    if (size(data,5)>1)
        if (size(data,5)~=dimensions.NoEchoes)
            return;
        end
    end
    if (size(data,6)>1)
        if (size(data,6)~=dimensions.NoExperiments)
            return;
        end
    end
    
    header1 = zeros(128,1); % 4x128=512 bytes

    header1(1)  = dimensions.NoSamples;
    header1(2)  = dimensions.NoViews;
    header1(3)  = dimensions.NoViews2;
    header1(4)  = dimensions.NoSlices;
    header1(39) = dimensions.NoEchoes;
    header1(40) = dimensions.NoExperiments;

    % Set datatype - 'complex float'
    header1(5)  = hex2dec('150000');

    % Open new file for writing
    fid1 = fopen(filename,'wb');

    % Write 512 byte header
    fwrite(fid1,header1,'int32');
    
    % For 3D data flip the 2nd and 3rd dimension
    if (size(data,3)>1)
        data = flip(permute(data,[2,3,1,4,5,6]),3);
    else
        data = flip(permute(data,[2,1,3,4,5,6]),2);
    end

    % Convert to 1D array with alternating real and imag part of the data
    temp = data;
    temp = temp(:);
    a = real(temp);
    b = imag(temp);
    temp = transpose([a b]);
    temp = temp(:);

    % Write data at once
    fwrite(fid1,temp,'float32');

    % write the footer
    fwrite(fid1,footer,'int8');

    % close file 
    fclose(fid1);

end