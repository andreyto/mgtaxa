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


## Recursive function that makes a path to a given directory traversable by everyone.
## It tries to explore both paths through symlinks and real paths.
_make_readable() {
    local opt
    local dir
    local prev_loop_dir
    local _debug
    _debug=$1
    for opt in -L -P; do
        dir=$(pwd $opt)
        #comparing to the dir processed with -L cuts down something
        #on the number of recursive calls when there are no links;
        #there is still excessive scanning that can be optimized but
        #this probably does not make a difference for any practical use.
        if [[ "$dir" != "/" && "$dir" != "$prev_loop_dir" ]]; then
            pushd $dir > /dev/null
            [ -n "$_debug" ] && echo "$opt Before: $(ls -ld $dir)"
            chmod +X . > /dev/null &> /dev/null
            [ -n "$_debug" ] && echo "$opt After: $(ls -ld $dir)"
            pushd .. > /dev/null
            _make_readable $_debug
            popd > /dev/null
            popd > /dev/null
        fi
        prev_loop_dir="$dir"
    done
}

## Make the target world-readable and directories leading to it traversable.
## If ther eare symlinks in the lineage, make their real lineages traversable as well.
make_readable() {
    local _debug
    _debug=1
    targ=$1
    [ -n "$targ" ] || exit 1
    pushd $(dirname $targ) > /dev/null
    chmod -R +Xr $(basename $targ)
    [ -n "$_debug" ] && echo "Before: $(ls -ld $targ)"
    _make_readable $_debug
    [ -n "$_debug" ] && echo "After: $(ls -ld $targ)"
}

