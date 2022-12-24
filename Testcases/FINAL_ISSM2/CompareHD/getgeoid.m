function geoidname = getgeoid(param)

%FUNCTION geoidname = getgeoid.m
%
% NAME: 
% VERSION: 1.0
%    
% PURPOSE: get ellisoid name from parameters vector (major axis and
% excentricity)
%    
% CATEGORY: map routines
%    
% CALLING SEQUENCE:
%    
% INPUTS:     param = vector of ellispoid parameters [major axis,excentricity] (in m)
%         
% OUTPUTS:  geoidname
%
%
% AUTHOR: Pascal Sirguey, University of Otago
%
% COPYRIGHT, LICENSE, and WARRANTY : refer to the MODImLAB User's Guide.

if isequal(param, almanac('earth','everest','meter'))
    geoidname = 'everest';
elseif isequal(param, almanac('earth','bessel','meter'))
    geoidname = 'bessel';
elseif isequal(param, almanac('earth','airy','meter'))
    geoidname = 'airy';
elseif isequal(param, almanac('earth','clarke66','meter'))
    geoidname = 'clarke66';
elseif isequal(param, almanac('earth','clarke80','meter'))
    geoidname = 'clarke80';
elseif isequal(param, almanac('earth','iau65','meter'))
    geoidname = 'iau65';
elseif isequal(param, almanac('earth','wgs66','meter'))
    geoidname = 'wgs66';
elseif isequal(param, almanac('earth','iau68','meter'))
    geoidname = 'iau68';
elseif isequal(param, almanac('earth','wgs72','meter'))
    geoidname = 'wgs72';
elseif isequal(param, almanac('earth','grs80','meter'))
    geoidname = 'grs80';
elseif isequal(param, almanac('earth','wgs84','meter'))
    geoidname = 'wgs84';
end
