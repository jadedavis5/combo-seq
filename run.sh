#!/usr/bin/env bash

# =========================== #
# - Nextflow Wrapper Script - #
# =========================== #
#
# Author - Yutathkarn Coles
# Email  - danielcoles1605@gmail.com
#
# Wrapper script written in Bash that runs Nextflow with default parameters,
# including loading profile-specific configs and parameters
# 

# Load a Nextflow profile, valid options are local and setonix
PROFILE="local"

case $PROFILE in
	(local)
		CONF="conf/local.conf"
		PARAMS="conf/local.yaml"
		;;
	(setonix)
		CONF="conf/setonix.yaml"
		PARAMS="conf/setonix.yaml"
		;;
	esac
