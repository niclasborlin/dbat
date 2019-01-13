# This is the README file for the Damped Bundle Adjustment Toolbox
# v0.8.5.1 for Matlab (R), below called the Toolbox.
#
# == LICENSE (short version) ==
#
# The Toolbox is released under a variant of the BSD license. The full
# text of the license is found in LICENSE.txt. In summary, you use the
# code at your own risk and may use it for any purpose, including
# commercial, as long as you give due credit. Furthermore, if you use
# the code, or derivatives thereof, for scientific publications, you
# should refer to one or more papers listed in the REFERENCES section.
#
#
# == INSTALLATION ==
#
# You can either install DBAT by downloading the source code or (if
# you use a git client) by cloning the repository.
#
# === Download ===
#
# 1) Download the package file dbat-master.zip (from the main page) or
#    dbat-x.y.z.w.zip/dbat-x.y.z.w.tar.gz (from the releases page) of
#    https://github.com/niclasborlin/dbat/
#
# 2) Unpack the file into a directory, e.g. c:\dbat or ~/dbat.
#
# === Clone ===
#
# At the unix/windows command line, write:
#
#   git clone https://github.com/niclasborlin/dbat.git
#
# to clone the repository into the directory 'dbat'. Use
#
#   git clone https://github.com/niclasborlin/dbat.git <dir-name>
#
# to clone the repository to another directory.
#
# If you use a graphical git client, e.g., tortoisegit
# (https://tortoisegit.org), select Git Clone... and enter
# https://github.com/niclasborlin/dbat.git or
# git@github.com:niclasborlin/dbat.git as the URL.
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
# == TESTING THE INSTALLATION ==
#
# 1) Start Matlab. Inside Matlab, do the following initialization:
# 1.1) cd c:\dbat % (change to where you unpacked the files)
# 1.2) dbatSetup  % will set the necessary paths, etc.
#
# 2) To test the demos, do 'help dbatdemos' or consult the manual.
#
#
# == UPDATING THE INSTALLATION==
#
# === Git ===
#
# If you cloned the archive, updating to the latest release is a
# simple as (replace ~/dbat and c:\dbat with where you cloned the
# repository):
#
#   cd ~/dbat
#   git pull
#
# at the command line. In TortoiseGit, right-click on the folder
# c:\dbat, select Git Sync... followed by Pull.
#
# === Download ===
#
# If you downloaded the code, repeat the download process under
# INSTALLATION. Most of the time it should be ok to unzip the new
# version on top of the old. However, we suggest you unzip the new
# version into a new directory, e.g. dbat-x-y-z-w, where x-y-z-w is
# the version number.
#
#
# == NEWS ==
# 
# For a list of recent changes, consult the ChangeLog.txt file.
# Updates will also be posted on https://www.researchgate.net/project/DBAT-The-Damped-Bundle-Adjustment-Toolbox-in-Matlab.
#
#
# == USAGE ==
#
# For examples of usage, see the usage section of the manual found in
# the doc/manual directory.
#
#
# == BUG REPORTS and/or FEATURE REQUESTS ==
#
# Instruction on how to submit bug reports and feature requests can be
# found in the file BUGREPORTS.txt
#
#
# == LICENSE ==
#
# Copyright (C) 2013-2019, Niclas Börlin, niclas.borlin@cs.umu.se (*),
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
#     DOI: 10.1111/phor.12037.
# Börlin and Grussenmeyer (2013). "Experiments with Metadata-derived
#     Initial Values and Linesearch Bundle Adjustment in Architectural
#     Photogrammetry", ISPRS Annals of the Photogrammetry, Remote
#     Sensing, and Spatial Information Sciences, vol II-5/W1:43-48.
#     DOI: 10.5194/isprsannals-II-5-W1-43-2013.
# Börlin and Grussenmeyer (2014). "Camera Calibration using the Damped
#     Bundle Adjustment Toolbox", ISPRS Annals of the Photogrammetry,
#     Remote Sensing, and Spatial Information Sciences, vol
#     II-5:89-96. DOI: 10.5194/isprsannals-II-5-89-2014.
# Börlin and Grussenmeyer (2016). "External Verification of the Bundle
#     Adjustment in Photogrammetric Software using the Damped Bundle
#     Adjustment Toolbox", ISPRS Archives XLI-B5, p. 7-14.
#     DOI: 10.5194/isprs-archives-XLI-B5-7-2016.
# Murtiyoso, Grussenmeyer, and Börlin (2017). "Reprocessing close
#     range terrestrial and UAV photogrammetric projects with the DBAT
#     toolbox for independent verification and quality control.",
#     ISPRS Archives XLII-2/W8:171-177.
#     DOI: 10.5194/isprs-archives-XLII-2-W8-171-2017.
# Börlin, Murtiyoso, Grussenmeyer, Menna, and Nocerino (2018).
#     "Modular bundle adjustment for photogrammeric computations",
#     ISPRS Archives XLII-2:133-140.
#     DOI: 10.5194/isprs-archives-XLII-2-133-2018.
#
#
# == CONTRIBUTIONS ==
#
# The toolbox includes contributions from the MathWorks File Exchange
# and the following people:
# * Arnaud Durand, ICube-SERTIT, University of Strasbourg, France.
# * Jan Hieronymus, TU Berlin, Germany.
# * Jean-Francois Hullo, EDF, France.
# * Fabio Menna, Fondazione Bruno Kessler, Trento, Italy.
# * Arnadi Murtiyoso, ICube, INSA Strasbourg, France.
# * Kostas Naskou, University of Nottingham, UK.
# * Erica Nocerino, Fondazione Bruno Kessler, Trento, Italy.
# * Deni Suwardhi, Bandung Institute of Technology, Indonesia.
#
#
# == TRADEMARKS ==
#
# Matlab is a registered trademark by The Mathworks, Inc., Natick MA,
# USA. PhotoModeler is a trademark owned by Eos System, Inc.,
# Vancouver, Canada. PhotoScan and Metashape are trademarks owned by
# Agisoft LLC, St. Petersburg, Russia. Other trademarks are owned by
# their respective owners.
