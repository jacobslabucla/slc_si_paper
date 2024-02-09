#!/bin/bash

# Path to the TSV file containing the list of directories
TSV_FILE="SMT_Neg_Files.tsv"

# Path to the new folder where you want to move the directories
NEW_FOLDER="../quant_SMT_Neg"

# Create the new folder if it doesn't exist
mkdir -p "$NEW_FOLDER"

# Read the TSV file line by line
while IFS=$'\t' read -r directory_name; do
# Print the directory name being processed
    echo "Processing directory: $directory_name"
   

    # Check if the directory exists
    if [ -d "$directory_name" ]; then
        # Move the directory to the new folder
        cp -r "$directory_name" "$NEW_FOLDER/"
        echo "Moved directory: $directory_name"
    else
        echo "Directory not found: $directory_name"
    fi
done < "$TSV_FILE"
