#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=3:00:00
#SBATCH --mail-type=FAIL
#SBATCH --job-name liam-omp-weak

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
args="-i $FINAL/data/squares.jpg -T -s 5 -M 2"
#Test epsilon starting from 6

function run_program() {
    cores=$1
    wk=$2
    LOGFILE="$FINAL/out/omp_weak_${cores}cores_eps${wk}.log"

    set -x
    # Run the program and save full output to log file
    # Extract time from output
    $EXE $args -c $cores -e $2 2>&1 | tee "$LOGFILE" | grep "tock" | cut -d' ' -f3 | tail -n1
    
    # Move output image to out directory if it exists
    if [ -f "$FINAL/img.jpg" ]; then
        mv "$FINAL/img.jpg" "$FINAL/out/omp_weak_${cores}cores_eps${wk}.jpg" 2>/dev/null || true
    fi
}

function run_with_num_cores() {
    echo "|core|eps|time|" 
    maxcores=40
    for i in "$@"; do
        wk=$(echo "scale=10; 6 - ((6 / $maxcores) * ($i - 1))" | bc)
        t1=$(run_program $i $wk)
        echo "|$i|$wk|$t1|"
    done
}

run_with_num_cores $(seq 1 40)
