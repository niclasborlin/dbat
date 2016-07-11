# Copy the installation section to a separate file.
perl -ne 'if (/^# == INSTALLATION/../^# == [^I]/) { print unless /^# == [^I]/; }' README.txt > ! INSTALL.txt
# Copy the license section to a separate file.
perl -ne 'if (/^# == LICENSE ==/../^# == [^LR]/) { print unless /^# == [^LR]/; }' README.txt > ! LICENSE.txt
