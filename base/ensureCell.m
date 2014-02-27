function r = ensureCell(z)
    %ensureCell Makes z into a cellarray if it is not one already. 
    %  Allows functions to easily accept either a single element or a cellarray of
    %  those elements as a parameter.
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    if(iscell(z))
        r = z;
    else
        r = {z};
    end
end
