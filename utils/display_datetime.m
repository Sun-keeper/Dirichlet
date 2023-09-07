function datetime_str = display_datetime(datetime_posix)
%DISPLAY_DATETIME Convert dates and times from posix format to human-readable formatted strings
%   datetime_str = DISPLAY_DATETIME(datetime_posix)
%
%   Input: a single numeric value of date and time in posix format or an array of these values
%
%   Output: a single string of date and time in the human-readable format or an array of these strings
%
% Author: Nhan Dam - March 30, 2018

datetime_obj = datetime(datetime_posix, 'ConvertFrom', 'posixtime');
datetime_str = datestr(datetime_obj);

end

