#!/bin/csh -f

#------------------------------------------------------------------------------------
# This script generates the antitracer_tavg_contents file for POP.
#
# It reads DIC and ALK master indices from antitracer_indices_dic.txt and
# antitracer_indices_alk.txt and generates the appropriate tavg variable names
# (e.g., DELTADIC007, DELTAALK007).
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
  echo "Error: Invalid my_stream number ($my_stream)"
  exit 5
endif

@ s1 = 1  # Use base-model stream 1

set dic_indices_file = "$CASEROOT/antitracer_indices_dic.txt"
set alk_indices_file = "$CASEROOT/antitracer_indices_alk.txt"
set output_file = "$CASEROOT/Buildconf/popconf/antitracer_tavg_contents"

if (! -e ${dic_indices_file}) then
    echo "Error: Required index file not found at ${dic_indices_file}"
    exit 7
endif

if (! -e ${alk_indices_file}) then
    echo "Error: Required index file not found at ${alk_indices_file}"
    exit 7
endif

set dic_indices = `cat ${dic_indices_file}`
set alk_indices = `cat ${alk_indices_file}`

# Create the file and write the tracer-independent diagnostics
cat >! ${output_file} << EOF
# Common antitracer-related gas exchange diagnostics
$s1  ANTITRACER_IFRAC
$s1  ANTITRACER_XKW
$s1  ANTITRACER_SCHMIDT
$s1  ANTITRACER_PV
EOF

if ($#dic_indices > 0) then
  foreach index ($dic_indices)
    set padded_index = `printf "%03d" $index`
    echo "$s1  DELTADIC${padded_index}"          >> ${output_file}
    echo "$s1  DELTADIC${padded_index}_FORCING"  >> ${output_file}
    echo "$s1  STF_DELTADIC${padded_index}"      >> ${output_file}
  end
endif

if ($#alk_indices > 0) then
  foreach index ($alk_indices)
    set padded_index = `printf "%03d" $index`
    echo "$s1  DELTAALK${padded_index}"          >> ${output_file}
    echo "$s1  DELTAALK${padded_index}_FORCING"  >> ${output_file}
    echo "$s1  STF_DELTAALK${padded_index}"      >> ${output_file}
  end
endif
