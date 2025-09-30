#!/bin/csh -f

#------------------------------------------------------------------------------------
# This script generates the antitracer_tavg_contents file for POP.
#
# It reads a list of master indices from antitracer_indices.txt and uses
# them to generate the specific tavg variable names (e.g., ANTITRACER007).
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
  echo "Error: Invalid my_stream number ($my_stream)"
  exit 5
endif

@ s1 = 1  # Use base-model stream 1

set indices_file = "$CASEROOT/antitracer_indices.txt"
set output_file = "$CASEROOT/Buildconf/popconf/antitracer_tavg_contents"

# Check that the file containing the master indices exists
if (! -e ${indices_file}) then
    echo "Error: Required index file not found at ${indices_file}"
    exit 7
endif

# Read the space-separated master indices from the file into an array
set master_indices = `cat ${indices_file}`

# Create the file and write the tracer-independent diagnostics
cat >! ${output_file} << EOF
# Common antitracer-related gas exchange diagnostics
$s1  ANTITRACER_IFRAC
$s1  ANTITRACER_XKW
$s1  ANTITRACER_SCHMIDT
$s1  ANTITRACER_PV
EOF

# Loop through the master indices and append the tracer-specific variables
foreach index ($master_indices)
  # Format the index to have three digits with leading zeros (e.g., 7 -> 007)
  set padded_index = `printf "%03d" $index`

  echo "$s1  ANTITRACER${padded_index}"           >> ${output_file}
  echo "$s1  ANTITRACER${padded_index}_FORCING"  >> ${output_file}
  echo "$s1  STF_ANTITRACER${padded_index}"      >> ${output_file}
  echo "$s1  ANTITRACER${padded_index}_COL_INT"  >> ${output_file}
end
