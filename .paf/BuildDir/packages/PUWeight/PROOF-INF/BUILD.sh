#!/bin/sh
# Build libEvent library.

VERBOSE=0

if [ "$VERBOSE" = "1" ]; then
    echo "#################### PUWeight #################################"
    echo "I am at `hostname` in dir `pwd`"
fi
 
if [ "$1" = "clean" ]; then
    [ "$VERBOSE" = "1" ] && echo "Cleaning..."
    [ "$VERBOSE" = "1" ] && make clean || make clean > /dev/null
    exit 0
fi

[ "$VERBOSE" = "1" ] && make || make > /dev/null
rc=$?

[ "$VERBOSE" = "1" ] && echo "rc=$?"
exit $rc

