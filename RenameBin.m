%script to rename .bin files resulting from microscope measurements
%be careful not to leave any other file inside the selected folder!!! 

path(path,'C:\Users\marti\OneDrive - Politecnico di Milano\MAGISTRALE\Thesis\Measurements');
dir0=('C:\Users\marti\OneDrive - Politecnico di Milano\MAGISTRALE\Thesis\Measurements');

dir2=uigetdir(dir0,'Select Directory');

el=dir(dir2);
fine=size(el,1)-2; % number of files in folder

ii=1;
while (el(ii).isdir)
    ii=ii+1;
end;

lun=length(el(ii).name); % Length of the shortest file

pnt=0;

h=waitbar(0,'Loading all files');

for ii=1:size(el,1)
    
    if not(el(ii).isdir)
        
        pnt=pnt+1;
        
        if length(el(ii).name)==lun
            
            oldname=[dir2,'/',el(ii).name];
            newname=[dir2,'/',el(ii).name(1:end-5),'000',el(ii).name(end-4:end)];
            % newname=[dir2,'/',el(ii).name(1:end-1),'00',el(ii).name(end:end)]; % ONLY for files with NO extension
            movefile (oldname,newname); % (1)
            
        elseif length(el(ii).name)==lun+1
            
            oldname=[dir2,'/',el(ii).name];
            newname=[dir2,'/',el(ii).name(1:end-6),'00',el(ii).name(end-5:end)];
            % newname=[dir2,'/',el(ii).name(1:end-2),'0',el(ii).name(end-1:end)]; % ONLY for files with NO extension
            movefile (oldname,newname); % (1)
            
        elseif length(el(ii).name)==lun+2
            
            oldname=[dir2,'/',el(ii).name];
            newname=[dir2,'/',el(ii).name(1:end-7),'0',el(ii).name(end-6:end)];
            % newname=[dir2,'/',el(ii).name(1:end-2),'0',el(ii).name(end-1:end)]; % ONLY for files with NO extension
            movefile (oldname,newname); % (1)
            
        end;
        
        waitbar(pnt/fine,h);
        
    end;
    
end

close(h);