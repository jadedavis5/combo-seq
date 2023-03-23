#!/usr/bin/env bash

# ===========================
# - Nextflow Wrapper Script -
# ===========================
#
# Author - Yutathkarn Coles
# Email  - danielcoles1605@gmail.com
#
# Wrapper script written in Bash that runs Nextflow with default parameters,
# including loading profile-specific configs and parameters
#

SCRIPT="workflow.nf"
NF_OPTS=""

# Python settings
# Set to a Python virtual environment, or comment out to disable
PYENV="containers/pipeline"
PY="python3"

# Enable auto-generation of Nextflow configs
AUTOGEN="True"
CONF_GEN="bin/generate_configs.py"
CONF_GET="bin/toml_parser.py"

# *** DO NOT EDIT BELOW THIS LINE ***

# Functions
function settings_get() {
	"$PY" "$CONF_GET" "$SETTINGS" "$1"
}

function test_python() {
	"$PY" --version &>/dev/null || echo "$NO_PYTHON"
	return 1
}

function test_toml() {
	"$PY" -c "import toml" &>/dev/null || echo "$NO_TOML"
	return 1
}

function load_venv() {
	echo "Loading Python virtual environment at $PYENV"
	source "$_PYENV" || echo "$NO_VENV"
}

function autogen() {
	# Relies on a functional Python 3 installation, and python3-toml
	# Python binaries
	CONF_DIR=$(settings_get "$_CONF")

	# Get location of Nextflow config and parameters
	CONFIG=$(settings_get "$_CONFIG")
	PARAMS=$(settings_get "$_PARAMS")

	# Set Nextflow settings
	if [[ $(settings_get "$_CLEAN_LOGS") == "$_BOOL_TRUE" ]]; then
		rm -rf .nextflow.log*
	fi
	if [[ $(settings_get "$_RESUME") == "$_BOOL_TRUE" ]]; then
		NF_OPTS+=" -resume"
	fi

	# Profile generation
	# Set active profile
	PROFILE=$(settings_get "$_PROFILE")

	# Generate configuration files from templates
	mkdir -p "$(dirname "$CONFIG")"
	mkdir -p "$(dirname "$PARAMS")"
	"$PY" "$CONF_GEN" "$CONF_DIR/$PROFILE.toml"
}

# Settings
# Run script paramters
SETTINGS="settings.toml"

# Settings strings
_PROFILE="profile.active"
_CONF="files.conf"
_CONFIG="files.config"
_PARAMS="files.params"
_CLEAN_LOGS="nextflow.clean_logs"
_RESUME="nextflow.resume"
_PYENV="$PYENV/bin/activate"
_NF="files.nf"
_PY="files.py3"

# Data types
# Boolean string for true
_BOOL_TRUE="True"
_BOOL_FALSE="False"

# Error messages
NO_GEN="""
Unable to autogenerate pipeline configuration
This will likely result in an unsable pipeline, please ensure that Python 3
and python3-toml are installed, manually editing auto-generated configurations
is highly discouraged"""
NO_VENV="""
Unable to load Python 3 virtual environment, using host Python 3 installation
> Tried to run 'source $_PYENV'
"""
NO_PYTHON="""
Python 3 appears to be unusable
> Unable to run binary: $PY
$NO_GEN"""
NO_TOML="""
Python 3 appears to be unable to load the module 'toml'
$NO_GEN"""

# Binaries
# NF="$(settings_get "$_NF")"
# PY="$(settings_get $_PY)"

NF="bin/nextflow"

# Load Python 3 virtual environment, else use host Python 3
if [[ -n "$PYENV" ]]; then load_venv; fi

# Environment checks
# Check that Python 3 and python3-toml are functioning
if [[ $(test_python) -gt 0 ]]; then AUTOGEN="$_BOOL_FALSE"; fi
if [[ $(test_toml) -gt 0 ]]; then AUTOGEN="$_BOOL_FALSE"; fi

# Ensure that AUTOGEN was set correctly, using the longest case conversion, ever
if [[ "$(echo $AUTOGEN | tr '[:upper:]' '[:lower:]')" == "$(echo $_BOOL_TRUE | tr '[:upper:]' '[:lower:]')" ]]; then
	AUTOGEN="$_BOOL_TRUE"
fi
if [[ "$AUTOGEN" == "$_BOOL_TRUE" ]]; then
	echo "Autogenerating Nextflow configuration files"
	autogen
else
	# Load fallback environment
	CONFIG="conf/generated/nextflow.conf"
	PARAMS="conf/generated/params.yaml"
fi

# Execute pipeline
echo """
> Running Nextflow ComboSeq pipeline

- Using Nextflow configuration: $CONFIG
- Loading paramters file: $PARAMS
"""

"$NF" -c "$CONFIG" run "$SCRIPT" -params-file "$PARAMS" $NF_OPTS
