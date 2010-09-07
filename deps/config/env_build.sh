
#!/bin/bash

## Source this file to set up custom build environment at JCVI

if [ -z "$AT_ENV_BUILD_DONE" ]; then

if [ -z "$AT_ENV_RUN_DONE" ]; then
    echo "Environment error: You need to source first the .env_run.sh"
    exit 1
fi

# If python has not been built yet, do not define related vars
if [ -f "$PYTHON" ]; then

#Some variables that help with Python extension module building for
#those installers that do not use distutils directly
export PY_INC_DIR=`$PYTHON -c 'from distutils.sysconfig import *; print get_python_inc()'`
#get_config_var() is "build dependent"
export PY_LIB_DIR=`$PYTHON -c 'from distutils.sysconfig import *; print get_config_var("LIBPL")'`
export PY_LIB=`$PYTHON -c 'from distutils.sysconfig import *; print get_config_var("LIBRARY")'`

fi #$PYTHON

#For C and C++ builds
export CFLAGS="-I${INSTMACH}/include -I${INST}/include ${CFLAGS}"
export CXXFLAGS="-I${INSTMACH}/include -I${INST}/include ${CXXFLAGS}"
export CPPFLAGS="${CXXFLAGS}"

#The -R will tell the linker to hard-wire the paths to shared libs as seen at linking.
#But LD_LIBRARY_PATH is still the only way to use some packages with older build procedures
export LDFLAGS="-L$INST_LIB_MACH -Wl,-R$INST_LIB_MACH -L${INST_LIB} -Wl,-R${INST_LIB} ${LDFLAGS}"


export AT_ENV_BUILD_DONE=1

fi # [ -z "$AT_ENV_BUILD_DONE" ]

