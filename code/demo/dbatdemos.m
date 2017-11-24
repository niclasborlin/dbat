%DBATDEMOS List of DBAT demos.
%
% loadplotdemo - Load and plot a camera/OP network.
%    Call:
%        LOADPLOTDEMO('ROMA') or LOADPLOTDEMO('CAM')
%    Data files:
%        DBATHOME/code/demos/data/dbat
%
%
% romabundledemo - Run bundle on Roma data set.
%    Call:
%        ROMABUNDLEDEMO
%    Data files:
%        DBATHOME/code/demos/data/dbat
%
%
% camcaldemo - Camera calibration example.
%    Call:
%        CAMCALDEMO
%    Data files:
%        DBATHOME/code/demos/data/dbat
%
%
% Prague 2016:
%    Call:
%        PRAGUE2016_PM(L,ORIENT), where L is 'C1'-'C2' or
%            'S1'-'S4', and ORIENT is TRUE or FALSE, or
%        PRAGUE2016_PS(L), where L is 'S5'
%    Data files:
%        DBATHOME/code/demos/data/prague2016/cam (C1-C2)
%        DBATHOME/code/demos/data/prague2016/sxb (S1-S5)
%
% Photoscan post-processing:
%    Call:
%        PS_POSTPROC(PSZNAME,SLOCAL) to post-process a .psz
%        project. Use SLOCAL=TRUE (recommended) to process in
%        semi-local coordinates (translation, scaling, no
%        rotation). Use PSZNAME='' to use the Prague 2016 SXB data set.
%    Data files:
%        DBATHOME/code/demos/data/prague2016/sxb/psprojects
%
% stpierrebundledemo_ps - Run bundle on StPierre data set from Photoscan.
%    Call:
%        STPIERREBUNDLEDEMO_PS
%    Data files:
%        DBATHOME/code/demos/data/hamburg2017/stpierre
%
%

% Run help if function is called.
help(mfilename)
