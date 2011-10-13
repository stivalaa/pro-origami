#!/bin/sh

# create diff to patch TableauCreate with
#
# cd TableauCreate 
# patch -p6 < patchfile
#
diff --recursive -c --new-file --exclude-from=excludefile.TableauCreate ~/TableauCreate TableauCreate
