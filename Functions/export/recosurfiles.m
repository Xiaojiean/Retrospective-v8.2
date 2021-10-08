function recosurfiles(surpath,suffix,mrdfilename,rprfilename)


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
    

end