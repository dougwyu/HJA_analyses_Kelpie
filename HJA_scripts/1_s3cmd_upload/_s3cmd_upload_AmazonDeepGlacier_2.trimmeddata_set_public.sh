#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# set a file public on Amazon S3
# this was run on macOS
#######################################################################################
#######################################################################################

s3cmd ls s3://amazon-oregon-douglasyu
# s3://amazon-oregon-douglasyu/2019Sep_shotgun_2.trimmeddata.tar
s3cmd setacl --acl-public s3://amazon-oregon-douglasyu/2019Sep_shotgun_2.trimmeddata.tar

# message back from Amazon:
	# "s3://amazon-oregon-douglasyu/2019Sep_shotgun_2.trimmeddata.tar: ACL set to Public  [1 of 1]"

# the public URL then follows this format:
	# https://<bucket-name>.s3.amazonaws.com/<object or key name>
https://amazon-oregon-douglasyu.s3.amazonaws.com/2019Sep_shotgun_2.trimmeddata.tar
# or from Transmit
https://amazon-oregon-douglasyu.s3.eu-west-1.amazonaws.com/2019Sep_shotgun_2.trimmeddata.tar


s3cmd ls --list-md5 s3://amazon-oregon-douglasyu/
