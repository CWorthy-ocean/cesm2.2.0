#!/bin/csh -f

#------------------------------------------------------------------------------------
# For now, set streams manually. You must only set as many streams as are declared
#  in the tavg_nml section. For example, if there are three streams:
#  @ s1 = $my_stream
#  @ s2 = $s1 + 1
#  @ s3 = $s2 + 1
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
   echo invalid my_stream number  ($my_stream)
   exit 5
endif

@ s1 = 1   # use base-model stream 1

cat >! $CASEROOT/Buildconf/popconf/antitracer_tavg_contents << EOF
$s1  ANTITRACER_IFRAC
$s1  ANTITRACER_XKW
$s1  ANTITRACER_SCHMIDT
$s1  ANTITRACER_PV
EOF

endif
