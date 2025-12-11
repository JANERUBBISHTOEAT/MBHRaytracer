#!/bin/bash
# Image comparison and difference map generation script
# Usage: ./compare_images.sh <image1> <image2> [output_prefix]

set -e

if [ $# -lt 2 ]; then
    echo "Usage: $0 <image1> <image2> [output_prefix]"
    echo "Example: $0 img.jpg data/squares_raytraced.jpg"
    exit 1
fi

IMAGE1="$1"
IMAGE2="$2"
OUTPUT_PREFIX="${3:-diff_map}"

if [ ! -f "$IMAGE1" ]; then
    echo "Error: $IMAGE1 not found"
    exit 1
fi

if [ ! -f "$IMAGE2" ]; then
    echo "Error: $IMAGE2 not found"
    exit 1
fi

# Get dimensions of second image for resizing
IMG2_DIMS=$(identify -format "%wx%h" "$IMAGE2")
IMG1_DIMS=$(identify -format "%wx%h" "$IMAGE1")

echo "Image 1: $IMAGE1 ($IMG1_DIMS)"
echo "Image 2: $IMAGE2 ($IMG2_DIMS)"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

# Create temporary directory
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Resize images to same size
echo "Resizing images..."
convert "$IMAGE1" -resize "$IMG2_DIMS" "$TMPDIR/img1_resized.jpg"
convert "$IMAGE2" "$TMPDIR/img2_resized.jpg"

# Calculate metrics
echo "Calculating comparison metrics..."
PSNR=$(compare -metric PSNR "$TMPDIR/img1_resized.jpg" "$TMPDIR/img2_resized.jpg" /dev/null 2>&1 || true)
RMSE=$(compare -metric RMSE "$TMPDIR/img1_resized.jpg" "$TMPDIR/img2_resized.jpg" /dev/null 2>&1 || true)
AE=$(compare -metric AE "$TMPDIR/img1_resized.jpg" "$TMPDIR/img2_resized.jpg" /dev/null 2>&1 || true)

# Generate difference maps
echo "Generating difference maps..."

# Grayscale difference map
convert "$TMPDIR/img1_resized.jpg" "$TMPDIR/img2_resized.jpg" \
    -compose difference -composite \
    -normalize -colorspace Gray \
    "${OUTPUT_PREFIX}.jpg"

# Color difference map
convert "$TMPDIR/img1_resized.jpg" "$TMPDIR/img2_resized.jpg" \
    -compose difference -composite \
    -normalize \
    "${OUTPUT_PREFIX}_color.jpg"

# Grayscale difference map (duplicate for consistency)
convert "$TMPDIR/img1_resized.jpg" "$TMPDIR/img2_resized.jpg" \
    -compose difference -composite \
    -normalize -colorspace Gray \
    "${OUTPUT_PREFIX}_gray.jpg"

# Generate side-by-side comparison
echo "Generating side-by-side comparison..."
convert \( "$TMPDIR/img1_resized.jpg" -label "$(basename "$IMAGE1") (resized)" \) \
    \( "$TMPDIR/img2_resized.jpg" -label "$(basename "$IMAGE2")" \) \
    \( "${OUTPUT_PREFIX}.jpg" -label "Difference Map" \) \
    +append -border 10x10 -bordercolor white \
    "${OUTPUT_PREFIX}_comparison.jpg"

# Calculate pixel-level statistics
echo "Calculating pixel-level statistics..."
STATS=$(convert "$TMPDIR/img1_resized.jpg" "$TMPDIR/img2_resized.jpg" \
    -compose difference -composite -colorspace Gray \
    -format "Total pixels: %[fx:w*h]\nMin difference: %[fx:minima*255]\nMax difference: %[fx:maxima*255]\nMean difference: %[fx:mean*255]\nStd deviation: %[fx:standard_deviation*255]\nMean difference %%: %[fx:mean*100]\nMax difference %%: %[fx:maxima*100]\nIdentical pixels: %[fx:mean==0?w*h:0] (%[fx:mean==0?100:0]%%)\n" \
    info: 2>&1)

# Channel-wise analysis
echo "Performing channel-wise analysis..."
convert "$TMPDIR/img1_resized.jpg" "$TMPDIR/img2_resized.jpg" \
    -compose difference -composite \
    -separate "$TMPDIR/diff_channels.png"

CHANNEL_STATS=""
for i in 0 1 2; do
    CHANNEL=$(echo "R G B" | awk -v i=$i '{print $((i+1))}')
    CH_STAT=$(convert "$TMPDIR/diff_channels-${i}.png" \
        -format "Channel $CHANNEL: mean=%[fx:mean*255], max=%[fx:maxima*255], stddev=%[fx:standard_deviation*255]\n" \
        info: 2>/dev/null)
    CHANNEL_STATS="${CHANNEL_STATS}${CH_STAT}"
done

# Output report
echo ""
echo "=== Image Comparison Report ==="
echo ""
echo "Files:"
echo "  Image 1: $IMAGE1 ($IMG1_DIMS)"
echo "  Image 2: $IMAGE2 ($IMG2_DIMS)"
echo ""
echo "Quality Metrics:"
echo "  PSNR: $PSNR dB"
echo "  RMSE: $RMSE"
echo "  Absolute Error Pixels: $AE"
echo ""
echo "Pixel-level Statistics:"
echo "$STATS"
echo ""
echo "Channel-wise Analysis:"
echo "$CHANNEL_STATS"
echo ""
echo "Generated Files:"
echo "  ${OUTPUT_PREFIX}.jpg - Grayscale difference map"
echo "  ${OUTPUT_PREFIX}_color.jpg - Color difference map"
echo "  ${OUTPUT_PREFIX}_gray.jpg - Grayscale difference map"
echo "  ${OUTPUT_PREFIX}_comparison.jpg - Side-by-side comparison"
echo ""
echo "Done!"
