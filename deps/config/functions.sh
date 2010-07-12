#!/bin/bash

# get_make_var VAR will find in the input stream VAR = VALUE and print VALUE
function get_make_var {
    sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*\(.*\)[[:space:]]*\(#.*\)*/\1/p"
}

# replace_make_var VAR VALUE will replace in the input stream the existing VAR = OLDVALUE with
# VAR = VALUE and print everything in the modified input stream
function replace_make_var {
    _par_var="$1"
    shift
    _par_val="$@"
    #escape / for sed into \/, using bash string subst
    _par_val=${_par_val//\//\\\/}
    sed "s/^\([[:space:]]*$_par_var[[:space:]]*=[[:space:]]*\)\(.*\)[[:space:]]*\(#.*\)*/\1$_par_val/"
}

