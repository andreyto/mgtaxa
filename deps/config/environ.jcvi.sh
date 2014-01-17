#!/bin/bash

## Source this file to set up custom environment needed e.g. by MGTAXA dependencies
## The intention is to keep the files .env_run.sh and .env_build.sh relatively generic,
## and put site-specific stuff in this file.

if [ -z "$MGT_PREFIX_ENVIRON_ENTERED" ]; then

## This will guard against cyclical sourcing in MGTAXA profile
export MGT_PREFIX_ENVIRON_ENTERED=1

this_file="$BASH_SOURCE"
## this file's abs dir
this_file_dir=$(cd $(dirname "$this_file") && pwd)

## this file's abs path
export MGT_PREFIX_ENVIRON_RC=$this_file_dir/$(basename "$this_file")

if [ -z "$MGT_PREFIX" ]; then
    ## Root of MGT installation - where all dependencies are
    ## In that case this file should be installed under packages/etc
    ## below the prefix
    export MGT_PREFIX=${this_file_dir%%/packages/etc}
fi

#export MGT_JCVI_DMZ=1

#Are we on system upgraded to CentOS 6?
if [[ $(uname -r) == *el6* ]]; then
    _cent6=1
else
    _cent6=
fi

## JCVI used to have the only has functional gcc (with gfortran) in /usr/local, so we need to make sure
## its shared libs are found at run-time instead of the one in /usr/
## Now it seems that gfortran has been installed from RPM in /usr, but all stuff in /usr/local probably
## was still compiled with gcc from /usr/local. Right now the OS and the /usr/local have the same
## version, so it still be OK if one shared lib is picked at run-time instead of another.
## On DMZ, the /usr/local gcc only lives under its own prefix, but we do not process DMZ case separately
## because the versions match with the OS gcc.
export CFLAGS="$CFLAGS -I/usr/local/include"
export LDFLAGS="$LDFLAGS -L/usr/local/lib64 -Wl,-R/usr/local/lib64 -L/usr/local/lib -Wl,-R/usr/local/lib"

if [ -z "$MGT_JCVI_DMZ" ]; then
# tell SGE to use new CentOS grid
source /usr/local/sge_current/jcvi/common/settings.sh
fi

export LD_LIBRARY_PATH=/usr/local/lib64:/usr/local/lib:$LD_LIBRARY_PATH

if [ -n "$_cent6" ]; then

    # Use conditional in case packages
    # have different locations between cluster and DMZ
    if [ -z "$MGT_JCVI_DMZ" ]; then
        true
    else
        true
    fi

    export JAVA_HOME=/usr/local/packages/jdk-6u41
    ## QPID needs java 6+
    export JAVA_6_HOME=$JAVA_HOME

    export GCC_HOME=/usr/local/packages/gcc-4.7.1

    #Now gcc is in the PATH but its libraries are not
    export PATH=$GCC_HOME/bin:$PATH
    export LD_LIBRARY_PATH=$GCC_HOME/lib64:$LD_LIBRARY_PATH
    export LDFLAGS="$LDFLAGS -L$GCC_HOME/lib64 -Wl,-R$GCC_HOME/lib64"

else

#This Java path is identical on Intranet and DMZ, although it seems that the main
#java home is under /usr/local/java on the Intranet
export JAVA_HOME=/usr/local/packages/java/current
## QPID needs java 6+
export JAVA_6_HOME=/usr/local/packages/java/1.6.0

fi

## ANT (the Java make) is installed at JCVI, but ANT_HOME has to be defined,
## otherwise Globus Toolkit build fails when ant cannot find xmlvalidate module ("task")
export ANT_HOME=/usr/local/packages/apache-ant

## We support these types of Python install:
## already built under a separate prefix and we do not have write support
## prefix same as the packages that we build, and we build our python where.
## We define the following variables to make it easier to support the above
## two cases and various other combinations.
## If set, we will build our own interpreter:
export PY_BUILD_INTERP=
## If set, the interpreter lives under that separate --prefix, otherwise
## it will be under one of common prefixes (either our package prefix or
## system prefix):
if [ -n "$_cent6" ]; then
    export PY_PREF_INTERP=/usr/local/packages/python-2.7.3
else
    export PY_PREF_INTERP=/usr/local/packages/python-2.6.6
fi
## Additionally, there is an option to set an absolute path to
## the Python executable (to support e.g. python2.5 and python2.6
## located inside the same bin/ directory). If not set, this variable
## will default to bin/python under the finally computed python prefix.
export PYTHON=

# now it is just an alias for MGT_PREFIX
# this has to point to the area with enough space, e.g. project area
export WORK="$MGT_PREFIX"

# Here we keep the source tree and the build area for all dependency packages
# Although this var is used only during build phase, it is more convenient to define
# it in this top level file.
# This will take 3G when the source tree is copied and then built.
# This can be deleted after everything has been built.
export DEP_SRC_TOP=$WORK/distros

# project area shared scratch space
if [ -z "$SCRATCH" ]; then
    if [ -z "$MGT_JCVI_DMZ" ]; then
        export SCRATCH=/usr/local/projects/MGTAXA/scratch/$USER
    else
        export SCRATCH=/opt/software/mgtaxa/scratch/$USER
    fi
fi

# This should be sourced at all times
source $MGT_PREFIX/packages/etc/env_run.sh
# This can be commented out after everything is built
source $MGT_PREFIX/packages/etc/env_build.sh

# Source MGTAXA config if exists
[ -f $INST/mgtaxa/etc/mgtaxa.shrc ] && source $INST/mgtaxa/etc/mgtaxa.shrc

# This gets rid of "Could not find/open font when opening font "arial", 
# using internal non-scalable font" in custom-compiled gnuplot
#export GDFONTPATH=/usr/share/fonts/liberation
#export GNUPLOT_DEFAULT_GDFONT=LiberationSans-Regular

export MGT_PREFIX_ENVIRON_DONE=1

fi # [ -z "$MGT_PREFIX_ENVIRON_DONE" ]


