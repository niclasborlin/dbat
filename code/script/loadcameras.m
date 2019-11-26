function [cam,xml]=loadcameras(fName)
%LOADCAMERAS Load camera specification from xml file.
%
%   CAM=LOADCAMERAS(CAMFILE) loads camera information from the DBAT
%   XML camera file CAMFILE. The camera information is returned in the
%   structure array CAM.
%
%See also: XMLTODBATCAMSTRUCT.

cam=[];

xml=dbatxml2struct(fName);

% Top level element is just 'document'.
[ok,msg]=checkxmlfields(xml,'document');

if ~ok
    error('DBAT camera XML top level error: %s',msg);
end

doc=xml.document;

% Document fields are input, output, operations, and Attributes.
docFields={'cameras','Attributes'};

[ok,msg]=checkxmlfields(doc,docFields);
if ~ok
    error('DBAT camera XML document field error: %s',msg);
end

% Known DBAT script version interval.
firstKnownVersion='1.0';
lastKnownVersion='1.0';

[ok,msg]=checkversionattr(doc.Attributes,'dbat_camera_version', ...
                          lastKnownVersion,firstKnownVersion);

if ~ok
    error('DBAT camera XML file error: %s',msg);
end

[ok,msg]=checkxmlfields(doc.cameras,'camera');
if ~ok
    error('DBAT camera XML document field error: %s',msg);
end

% Extract and return camera list.
cam=parsedbatxmlcamstruct(doc.cameras.camera);
