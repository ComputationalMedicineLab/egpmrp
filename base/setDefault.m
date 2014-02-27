function [] = setDefault( name, val )
    %SETDEFAULT Sets the value of variable named 'name' to val if name is empty or
    %non-existent.
    %
    % Intended to be used to initialize variables inside a function call if their
    % value was not passed to the function or was passed as an empty vector.
    % Usage:
    %
    % function val = f(x)
    %   setDefault('x', 3);
    %   % now ok to use x
    % end
    %
    % Note the quotes around x.
    %
    % Because MATLAB doesn't have lazy evaluation, usage such as
    % setDefault('x', ExpensiveInitializationFunction());
    % will  evaluate ExpensiveInitializationFunction() before calling setDefault,
    % regardless of whether x will be assigned its return value. To approximate lazy
    % evaluation and avoid the unnecessary computation or other side effects, you
    % can wrap the initialization in an anonymous function:
    %
    % setDefault('x', @() ExpensiveInitializationFunction());
    %
    % This will evaluate the function only if it is actually needed to initialize x.
    % If the function has been defined elsewhere, you can pass its handle,
    % which will also avoid the unnecessary computation:
    %
    % setDefault('x', @ExpensiveInitializationFunction);
    %
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    assert(ischar(name), 'Variable name must be a string.');
    
    if evalin('caller', ['exist(''', name, ''', ''var'') ~= 1 || isempty(', name, ')'])
        if isa(val, 'function_handle')
            assignin('caller', name, val());
        else
            assignin('caller', name, val);
        end
    end
    
end

