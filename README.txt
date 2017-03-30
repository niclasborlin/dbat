# This is the README file for the Damped Bundle Adjustment Toolbox
# v0.6.2.1 for Matlab (R), below called the Toolbox.
#
# == LICENSE (short version) ==
#
# The Toolbox is released under a variant of the BSD license. The full
# text of the license is found in LICENSE.txt. In summary, you use the
# code at your own risk and may use it for any purpose, including
# commercial, as long as you give due credit. Furthermore, if you use
# the code, or derivatives thereof, for scientific publications, you
# should refer to on or more papers listed in the REFERENCES section.
#
#
# == INSTALLATION ==
#
# 1) Download the package file dbat-x.y.z.zip or dbat-x.y.z.tar.gz
#    from https://github.com/niclasborlin/dbat/
#
# 2) Unpack the file into a directory, e.g. c:\dbat or ~/dbat.
# 
# 3) Start Matlab. Inside Matlab, do the following initialization:
# 3.1) cd c:\dbat % (change to where you unpacked the files)
# 3.2) dbatSetup  % will set the necessary paths, etc.
#
# 4) To test the demos, do 'help dbatdemos' or consult the manual.
#
#
# ==== Download high-resolution images ====
#
# To reduce the size of the repository and hence download times, only
# low-resolution images are included in the repository. High-resolution 
# images can be downloaded from http://www.cs.umu.se/~niclas/dbat_images/.
# For further details, consult the README.txt files in the respective
# image directories.
#
#
# == NEWS ==
# 
# For a list of recent changes, consult the ChangeLog.txt file.
#
#
# == USAGE ==
#
# For examples of usage, see the usage section of the manual found in
# the doc/manual directory.
#
#
# == LICENSE ==
#
# Copyright (C) 2013-2017, Niclas Börlin, niclas.borlin@cs.umu.se (*),
# and Pierre Grussenmeyer, pierre.grussenmeyer@insa-strasbourg.fr (**).
# All rights reserved.
#
# (*)  Department of Computing Science, Umeå University, Sweden. 
# (**) ICube Laboratory UMR 7357, Photogrammetry and Geomatics Group,
#      Institut national des sciences appliquées de Strasbourg,
#      Strasbourg, France.  
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#
# 3. If you use the code, or derivatives thereof, for scientific
#    publications, you should refer to on or more of the papers in the
#    REFERENCES section.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
#
# == REFERENCES ==
#
# Börlin and Grussenmeyer (2013). "Bundle adjustment with and without
#     damping", Photogrammetric Record, vol. 28(144):396-415.
# Börlin and Grussenmeyer (2013). "Experiments with Metadata-derived
#     Initial Values and Linesearch Bundle Adjustment in Architectural
#     Photogrammetry", ISPRS Annals of the Photogrammetry, Remote
#     Sensing, and Spatial Information Sciences, vol II-5/W1:43-48.
# Börlin and Grussenmeyer (2014). "Camera Calibration using the Damped
#     Bundle Adjustment Toolbox", ISPRS Annals of the Photogrammetry,
#     Remote Sensing, and Spatial Information Sciences, vol
#     II-5:89-96.
# Börlin and Grussenmeyer (2016). "External Verification of the Bundle
#     Adjustment in Photogrammetric Software using the Damped Bundle
#     Adjustment Toolbox", ISPRS Archives XLI-B5, p. 7-14.
#
#
# == CONTRIBUTIONS ==
#
# The toolbox includes contributions from the MathWorks File Exchange
# and the following people:
# * Arnaud Durand, ICube-SERTIT, University of Strasbourg.
# * Jan Hieronymus, TU Berlin.
# * Jean-Francois Hullo, EDF.
# * Kostas Naskou, University of Nottingham.
#
#
# == TRADEMARKS ==
#
# Matlab is a registered trademark by The Mathworks, Inc., Natick MA,
# USA.
#
