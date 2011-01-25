#!/bin/bash

## Source this file to set up custom run-time environment at JCVI

if [ -z "$AT_ENV_RUN_DONE" ]; then

# I would want to get rid of /usr/local/bin. Lots of broken stuff.
# But 'use localbins' somehow breaks sudo access. Sudo starts claiming that
# my account has no sudo rights (test with 'sudo -l')
#use localbins


## Some JCVI profile files apparently wack LD_LIBRARY_PATH
## Source MPI selector settings if it is installed on the local system
#[ -f /etc/profile.d/mpi-selector.sh ] && source /etc/profile.d/mpi-selector.sh

export INST=$WORK/packages

export CPUARCH=x86_64
export DISTRO=rhel5
export MACH=x86_64-rhel5

export INSTMACH=$INST/$MACH

export INST_LIB_MACH=$INSTMACH/lib64

export INST_LIB=$INSTMACH/lib

## The location of files with initialization data
## that should remain private and not go into a public repo 
## - such as initial PostgreSQL root password.
## This directory should be readable by user only
export PRIV_DATA=$HOME/.priv_data

if [ -n "$PY_PREF_INTERP" ]; then
## For python under its own prefix, we need to add paths to relevant env vars
export LD_LIBRARY_PATH="$PY_PREF_INTERP/lib:$LD_LIBRARY_PATH"
export PATH="$PY_PREF_INTERP/bin:$PATH"
[ -n "$PYTHON" ] || export PYTHON=$PY_PREF_INTERP/bin/python
else
## Otherwise, the common paths are already added
[ -n "$PYTHON" ] || export PYTHON=$INSTMACH/bin/python
fi

if [ -f "$PYTHON" ]; then

export PY_VER=`$PYTHON -c 'from distutils.sysconfig import *; print get_python_version()'`

export PYCOMMON=$INST/lib/python${PY_VER}/site-packages
export PYMACH=$INSTMACH/lib/python${PY_VER}/site-packages
# Shogun will install itself in here if given INSTMACH as prefix,
# and this is clumsy to override, so we just add it to the path:
export PYDIST=$INSTMACH/lib/python${PY_VER}/dist-packages
export PYTHONPATH=${PYDIST}:${PYMACH}:${PYCOMMON}:${PYTHONPATH}

fi

export PATH="$INST/x86_32/texlive/2007/bin/i386-linux:$INSTMACH/bin:$INST/bin:$PATH"


# For some packages, setting LD_RUN_TIME during build and linker options do not help.
# We still have to define LD_LIBRARY_PATH

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${INST_LIB_MACH}:${INST_LIB}:${INST_LIB_MACH}/R/lib:${INST_LIB}/mysql"

# This will be used when building our own MySQL package, and to configure MGTAXA connection
export MGT_MYSQL_HOST=mgtaxa-dev.jcvi.org
export MGT_MYSQL_PORT=13306
export MGT_MYSQL_TMPDIR=$SCRATCH

# Globus Toolkit
export GLOBUS_LOCATION=$INSTMACH/globus
export PATH=$GLOBUS_LOCATION/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GLOBUS_LOCATION/lib

# GridWay
export GW_LOCATION=$INSTMACH/gridway
export PATH=$GW_LOCATION/bin:$PATH
## We use compile-time linker flags instead:
#export LD_LIBRARY_PATH="${GW_LOCATION}/lib:${LD_LIBRARY_PATH}"

# GridWay DRMAA Python (also needs LD_LIBRARY_PATH to include 
export PYTHONPATH=$GW_LOCATION/lib/python:${PYTHONPATH}

# PostgreSQL
export PGSQL_LOCATION=$INSTMACH/pgsql
export PATH=$PGSQL_LOCATION/bin:$PATH
export PGSQL_DATA=$PGSQL_LOCATION/data
export PGDATA=$PGSQL_DATA # PostgreSQL recognizes this
export LD_LIBRARY_PATH="${PGSQL_LOCATION}/lib:${LD_LIBRARY_PATH}"
## PostgreSQL client psql will try to get passwords from this file -
## this is th way to avoid password prompts for maintenance scripts
export PGPASSFILE=$PRIV_DATA/pgpass
export PGSQL_ROOT_USER=postgres

# Apache Qpid
export QPID_VER=0.8
export QPID_HOME=$INSTMACH/qpid-$QPID_VER
export QPID_WORK=$QPID_HOME/var
export PATH=$QPID_HOME/bin:$PATH
export QPID_JAVA_HOME=/usr/local/packages/java/1.6.0

# Galaxy
export GALAXY_LOCATION=$INSTMACH/mgtaxa-galaxy

# Ruby stuff
# We need to have rubygems module loaded ("require") on every start of ruby,
# because that adds into Ruby's library path the gems installed by 'gem install --user-install'
# into user's home dir.
# Otherwise Ruby cannot find them.
# @todo We probably need to add GEM_HOME=$INSTMACH/ruby/gem so that --user-install puts gems
# where instead of default ~/.gem
export RUBYOPT="-rrubygems"

export AT_ENV_RUN_DONE=1

fi # [ -z "$AT_ENV_RUN_DONE" ]


