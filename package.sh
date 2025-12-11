#!/bin/bash
# Package the entire project directory into a zip file
# Usage: ./package.sh [output_filename]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_FILE="${1:-BHRaytracer.zip}"

# If output file has no path, put it in the script directory
if [[ "$OUTPUT_FILE" != */* ]]; then
    OUTPUT_FILE="$SCRIPT_DIR/$OUTPUT_FILE"
fi

echo "Packaging project directory..."
echo "Source: $SCRIPT_DIR"
echo "Output: $OUTPUT_FILE"
echo ""

# Remove old zip file if it exists
if [ -f "$OUTPUT_FILE" ]; then
    echo "Removing existing zip file..."
    rm -f "$OUTPUT_FILE"
fi

# Change to project directory
cd "$SCRIPT_DIR"

# Create zip file, excluding zip files themselves
echo "Creating zip archive..."
zip -r "$OUTPUT_FILE" . -x "*.zip" > /dev/null 2>&1

# Get file statistics
ZIP_SIZE=$(ls -lh "$OUTPUT_FILE" | awk '{print $5}')
ZIP_INFO=$(zipinfo "$OUTPUT_FILE" | tail -1)

echo "Package created successfully!"
echo ""
echo "File: $OUTPUT_FILE"
echo "Size: $ZIP_SIZE"
echo "Contents: $ZIP_INFO"
echo ""
