#!/bin/csh -f

#------------------------------------------------------------------------------------
# This script generates the antitracer_tavg_contents file for POP.
# It supports multiple antitracers based on the antitracer_tracer_cnt setting.
#
# Assumes 'ANTITRACER_TRACER_CNT' is available as an environment variable
# set by the CESM build system (e.g., in env_build.xml or through a CIME macro).
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
  echo "Error: Invalid my_stream number ($my_stream)"
  exit 5
endif

@ s1 = 1   # Use base-model stream 1 (or other relevant stream index)

# Get the number of antitracers.
# This variable (ANTITRACER_TRACER_CNT) MUST be set in your CESM environment
# (e.g., in env_build.xml or through a CIME macro during case setup).
# Example: setenv ANTITRACER_TRACER_CNT 2
if (! $?ANTITRACER_TRACER_CNT) then
    echo "Error: ANTITRACER_TRACER_CNT environment variable is not set."
    echo "Please ensure it's defined in your CESM case environment."
    exit 6
endif

@ num_antitracers = ${ANTITRACER_TRACER_CNT}

cat >! $CASEROOT/Buildconf/popconf/antitracer_tavg_contents << EOF
# Common antitracer-related gas exchange diagnostics
$s1  ANTITRACER_IFRAC
$s1  ANTITRACER_XKW
$s1  ANTITRACER_SCHMIDT
$s1  ANTITRACER_PV
EOF

# Loop through each antitracer to add its concentration and flux to the tavg_contents
@ n = 1
while ($n <= $num_antitracers)
    # Append the tracer's concentration
    echo "$s1  ${ANTITRACER_TRACER_NAMES[$n]}" >> $CASEROOT/Buildconf/popconf/antitracer_tavg_contents
    # Append the tracer's surface flux (assuming your Fortran module names it ANTITRACER_SFLUX or similar)
    # NOTE: You need to define the name for the surface flux if it's not simply the tracer name.
    # If the surface flux for Antitracer_N is named "ANTITRACER_N_SFLUX" in your Fortran output:
    echo "$s1  ${ANTITRACER_TRACER_NAMES[$n]}_SFLUX" >> $CASEROOT/Buildconf/popconf/antitracer_tavg_contents
    # If your tracer name is already `ANTITRACER1` and you want that as the variable, it works.
    # If you wanted a specific variable like the surface tracer itself or its flux:
    # Example for surface tracer value:
    # echo "$s1  TS_OCN_${ANTITRACER_TRACER_NAMES[$n]}" >> $CASEROOT/Buildconf/popconf/antitracer_tavg_contents

    @ n = $n + 1
end
