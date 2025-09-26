#!/bin/csh -f

#------------------------------------------------------------------------------------
# This script generates the antitracer_tavg_contents file for POP.
# It supports multiple antitracers based on the antitracer_tracer_cnt setting.
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
  echo "Error: Invalid my_stream number ($my_stream)"
  exit 5
endif

@ s1 = 1  # Use base-model stream 1

if (! $?ANTITRACER_TRACER_CNT) then
    echo "Error: ANTITRACER_TRACER_CNT environment variable is not set."
    exit 6
endif

@ num_antitracers = ${ANTITRACER_TRACER_CNT}

# Create the file and add the common diagnostics
cat >! $CASEROOT/Buildconf/popconf/antitracer_tavg_contents << EOF
# Common antitracer-related gas exchange diagnostics
$s1  ANTITRACER_IFRAC
$s1  ANTITRACER_XKW
$s1  ANTITRACER_SCHMIDT
$s1  ANTITRACER_PV
EOF

# Append the main tracer fields and their averages to the file
@ i = 1
while ($i <= $num_antitracers)
  echo "$s1  ANTITRACER${i}"          >> $CASEROOT/Buildconf/popconf/antitracer_tavg_contents
  echo "$s1  ANTITRACER${i}_FORCING" >> $CASEROOT/Buildconf/popconf/antitracer_tavg_contents
  echo "$s1  STF_ANTITRACER${i}"    >> $CASEROOT/Buildconf/popconf/antitracer_tavg_contents
  echo "$s1  ANTITRACER${i}_COL_INT"  >> $CASEROOT/Buildconf/popconf/antitracer_tavg_contents  # << ADD THIS LINE
  @ i++
end
