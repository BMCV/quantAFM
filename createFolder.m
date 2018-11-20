function [success] = createFolder(folderName)
%CREATEFOLDER Creates a folder if it doesn't already exist
%   Checks whether a given folder path exists, and - based on the result -
%   creates a named folder, or does nothing.
    if ~ exist(fullfile(folderName), 'dir')
            mkdir(folderName)
    end
    success = 1;
end

