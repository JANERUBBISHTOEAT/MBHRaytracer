#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=3:00:00
#SBATCH --mail-type=FAIL
#SBATCH --job-name liam-omp-strong

FINAL=$1
EXE=$1/openmp/main

source $FINAL/scripts/teachsetup

if [ ! -f $EXE ]; then
    echo "Forgot to compile final, was looking for: '$EXE'"
    echo "Maybe try adding an argument on your sbatch invocation"
    echo "sbatch path/to/script path/to/final"
    exit 1
fi

#Pretty grind-y settings.
# Use -M 2 for mass, other parameters use defaults
# Ensure dual black hole mode (no -1 parameter)
args="-i $FINAL/data/squares.jpg -T -s 3 -M 2"

function run_program() {
    cores=$1
    LOGFILE="$FINAL/out/omp_strong_${cores}cores.log"
    
    # Run the program and save full output to log file
    # Extract time from output
    $EXE $args -c $cores 2>&1 | tee "$LOGFILE" | grep "tock" | cut -d' ' -f3
    
    # Move output image to out directory if it exists
    if [ -f "$FINAL/img.jpg" ]; then
        mv "$FINAL/img.jpg" "$FINAL/out/omp_strong_${cores}cores.jpg" 2>/dev/null || true
    fi
}

function run_with_num_cores() {
    echo "|core|time|" 
    for i in "$@"; do
        t1=$(run_program $i)
        echo "|$i|$t1|"
    done
}

echo "Ready to begin"
run_with_num_cores $(seq 1 20)
echo "Done!"
