#!/bin/bash

# Total count of 0.00000000E+00 = NS-3
value=${1}
total_count=$((value-3))

file="tprbap_modules.f"
#file="test.txt"
echo "> Changing the value of NI in line 8 of $file to $total_count"
sed -i "8s/\(NI=\)[0-9]*/\1$total_count/" "$file"
echo "> Recompile after modifying $file"
make clean && make

# Remove existing fort.43 file
rm fort.43
file_name="fort.43"

# Number of items per line
per_line=5

# Calculate the number of lines
line_count=$((total_count / per_line))

# Calculate the remainder for the last line
remainder=$((total_count % per_line))

# Start generating the file
echo -n > "$file_name"

# Iterate to generate each line
for ((i = 0; i < line_count; i++)); do
    echo -n " " >> "$file_name"
    for ((j = 0; j < per_line; j++)); do
        echo -n " 0.00000000E+00 " >> "$file_name"
    done
    echo >> "$file_name"
done

# Generate the last line
echo -n " " >> "$file_name"
for ((i = 0; i < remainder; i++)); do
    echo -n " 0.00000000E+00 " >> "$file_name"
done

# Check if there's a remainder
if [[ $remainder != 0 ]]; then
    echo " " >> "$file_name"
    echo -n " " >> "$file_name"
    echo -n " " >> "$file_name"
else
    echo -n " " >> "$file_name"
fi

# Display the contents of the generated file
cat "$file_name"


