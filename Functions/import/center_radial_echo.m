function output_data = center_radial_echo(app,input_data)

% shifts the echoes to the center = dimx/2 + 1

[nr,dimy,dimx] = size(input_data);


factor = 8;
center = round(dimx*factor/2 + 1);


TextMessage(app,'Centering the radial echoes ...');

parfor i = 1:nr
    
    for j= 1:dimy
        
        data_int = interp(squeeze(input_data(i,j,:)),factor);
        
        m = find(abs(data_int) == max(abs(data_int)));
        shift = round((dimx*factor)/2 - m - factor/2);
        data_int = circshift(data_int,shift);
        
        input_data(i,j,:) = downsample(data_int,factor);
        
    end
    
end

TextMessage(app,'Phase correction ...');

data_int0 = interp(squeeze(input_data(1,1,:)),factor);
phi0 = angle(data_int0(center));

parfor i = 1:nr
    
    for j= 1:dimy
        
        data_int = interp(squeeze(input_data(i,j,:)),factor);
        phi1 = angle(data_int(center));
        dphi = phi1 - phi0;
        data_int = data_int*exp(-1i*dphi);
        
        input_data(i,j,:) = downsample(data_int,factor);
        
    end
    
end


output_data = input_data;

end

