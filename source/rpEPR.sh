#!/bin/bash
# A script to make fifos, launch gnuplot, and the rpEPR program
# This allows the rpEPR program to still work even when
# compiled without posix support ie. with gfortran

# Make the fifos and launch gnuplot
rm plotcmds my_mouse 2 &> 1 >/dev/null
mkfifo plotcmds 2 &> 1 >/dev/null
mkfifo my_mouse 2 &> 1 >/dev/null
gnuplot <plotcmds &

# The default buffered io with gfortran sits on the data
# and gnuplot doesn't see it until fortran closes the file!
export GFORTRAN_UNBUFFERED_ALL=y

# Launch rpEPR
/usr/local/bin/gnu_rpEPR "$@"

# Clean up
rm plotcmds 2 &> 1 >/dev/null
rm my_mouse 2 &> 1 >/dev/null
unset GFORTRAN_UNBUFFERED_ALL
