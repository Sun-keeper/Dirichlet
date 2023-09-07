function KeyPressFcnTest
%KEYPRESSFCNTEST Detect the code for a key press
%   How to use this function: call this function in Command Window. Then, a figure window appears. Now we can press any
%   key that we want to find its code and the results will be printed in Command Window. To stop this program, close the
%   figure window (by clicking the 'close' on the top corner of its window) or type in Command Window 'close all'/'close
%   <id>' (where <id> is the id of the figure window).
%
% Source: https://au.mathworks.com/matlabcentral/answers/247892-how-to-detect-esc-key-press

close all;
h = figure;

set(h,'WindowKeyPressFcn',@KeyPressFcn);

    function KeyPressFcn(~,evnt)
        fprintf('key event is: %s\n',evnt.Key);
    end

end

