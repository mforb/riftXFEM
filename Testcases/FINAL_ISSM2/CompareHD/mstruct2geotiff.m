function key = mstruct2geotiff(userprojstruct);
%FUNCTION 
%
% NAME: 
% VERSION: 1.2.5
%
% PURPOSE: 
%
% CATEGORY:MODImLab
%
% CALLING SEQUENCE:
%
% INPUTS:
%
% OUTPUTS:
%       key: GeoTIFFTags structure
%
% AUTHOR: Pascal Sirguey, University of Otago
%
% TRACK CHANGE from v1.2.5
%   v1.2.4: - initial release
%   v1.2.5: - implement stereo
%
% COPYRIGHT, LICENSE, and WARRANTY : refer to the MODImLAB User's Guide.

% 6.2.1 GeoTIFF Configuration Keys
   key.GTModelTypeGeoKey            = 1;     % ModelTypeProjected   
   key.GTRasterTypeGeoKey           = 1;     % RasterPixelIsArea  
   key.GTCitationGeoKey             = 'user defined';
   
%6.2.2 Geographic CS Parameter Keys
   key.GeographicTypeGeoKey         = 32767; % user-defined GCS
   key.GeogCitationGeoKey           = 'user defined';
   switch getgeoid(userprojstruct.geoid)
       case 'everest'
           key.GeogGeodeticDatumGeoKey      = 6015; % DatumE_Everest1830_1937Adjustment 
       case 'bessel'
           key.GeogGeodeticDatumGeoKey      = 6004; % DatumE_Bessel_1841 
       case 'airy'
           key.GeogGeodeticDatumGeoKey      = 6001; % DatumE_Airy1830 
       case 'clarke66'
           key.GeogGeodeticDatumGeoKey      = 6008; % DatumE_Clarke_1866
       case 'clarke80'
           key.GeogGeodeticDatumGeoKey      = 6011; % DatumE_Clarke_1880_IGN 
       case 'grs80'
           key.GeogGeodeticDatumGeoKey      = 6019; % DatumE_GRS_1980 
       case 'wgs84'
           key.GeogGeodeticDatumGeoKey      = 6030; % DatumE_WGS_84 
       otherwise
           key.GeogGeodeticDatumGeoKey      = 32767; % user-defined GCS           
   end
   key.GeogPrimeMeridianGeoKey      = 8901;  % PM_Greenwich 
   key.GeogLinearUnitsGeoKey        = 9001;  % Linear_Meter 
   %key.GeogLinearUnitSizeGeoKey     =
   key.GeogAngularUnitsGeoKey       = 9102;  % Angular_Degree 
   %key.GeogAngularUnitSizeGeoKey    = 
   switch getgeoid(userprojstruct.geoid)
       case 'everest'
           key.GeogEllipsoidGeoKey      = 7015; % Ellipse_Everest_1830_1937_Adjustment 
       case 'bessel'
           key.GeogEllipsoidGeoKey      = 7004; % Ellipse_Bessel_1841 
       case 'airy'
           key.GeogEllipsoidGeoKey      = 7001; % Ellipse_Airy_1830 
       case 'clarke66'
           key.GeogEllipsoidGeoKey      = 7008; % Ellipse_Clarke_1866
       case 'clarke80'
           key.GeogEllipsoidGeoKey      = 7011; % Ellipse_Clarke_1880_IGN 
       case 'grs80'
           key.GeogEllipsoidGeoKey      = 7019; % Ellipse_GRS_1980 
       case 'wgs84'
           key.GeogEllipsoidGeoKey      = 7030; % Ellipse_WGS_84 
       otherwise
           key.GeogEllipsoidGeoKey      = 32767; % user-defined GCS           
   end
   key.GeogSemiMajorAxisGeoKey      = userprojstruct.geoid(1);
   %key.GeogSemiMinorAxisGeoKey      = 
   key.GeogInvFlatteningGeoKey      = userprojstruct.geoid(1)/(userprojstruct.geoid(1)-sqrt(userprojstruct.geoid(1)^2*(1-userprojstruct.geoid(2)^2)));
   %key.GeogAzimuthUnitsGeoKey       = 
   key.GeogPrimeMeridianLongGeoKey  = 0;
   
% 6.2.3 Projected CS Parameter Keys     
   key.ProjectedCSTypeGeoKey          = 32767; % user-defined
   key.PCSCitationGeoKey              = 'user defined';
   key.ProjectionGeoKey               = 32767; % user-defined
   switch userprojstruct.mapprojection
       case 'tranmerc'
           key.ProjCoordTransGeoKey           = 1; % CT_TransverseMercator 
       case {'lambert','lambertstd'}
           key.ProjCoordTransGeoKey           = 8; % CT_LambertConfConic_2SP 
       case {'stereo'}
           key.ProjCoordTransGeoKey           = 15; % CT_PolarStereographic
   end
   key.ProjLinearUnitsGeoKey          = 9001;  % Linear_Meter 
   %key.ProjLinearUnitSizeGeoKey       = 
   switch userprojstruct.mapprojection
       case {'lambert','lambertstd'}
           key.ProjStdParallel1GeoKey         = userprojstruct.mapparallels(1);
           key.ProjStdParallel2GeoKey         = userprojstruct.mapparallels(2);
end
   
   key.ProjNatOriginLongGeoKey        = userprojstruct.origin(2);
   key.ProjNatOriginLatGeoKey         = userprojstruct.origin(1);
   key.ProjFalseEastingGeoKey         = userprojstruct.falseeasting;
   key.ProjFalseNorthingGeoKey        = userprojstruct.falsenorthing;
   %key.ProjFalseOriginLongGeoKey      = 
   %key.ProjFalseOriginLatGeoKey       = 
   %key.ProjFalseOriginEastingGeoKey   = 
   %key.ProjFalseOriginNorthingGeoKey  = 
   %key.ProjCenterLongGeoKey           = 
   %key.ProjCenterLatGeoKey            = 
   %key.ProjCenterEastingGeoKey        = 
   %key.ProjCenterNorthingGeoKey       = 
   key.ProjScaleAtNatOriginGeoKey     = userprojstruct.scalefactor;
   %key.ProjScaleAtCenterGeoKey        = 
   %key.ProjAzimuthAngleGeoKey         = 
   %key.ProjStraightVertPoleLongGeoKey = 
   
% 6.2.4 Vertical CS Keys   
   %key.VerticalCSTypeGeoKey           = 0;  % not defined
   %key.VerticalCitationGeoKey         = 'user defined';
   %key.VerticalDatumGeoKey            = 
   %key.VerticalUnitsGeoKey            = 
