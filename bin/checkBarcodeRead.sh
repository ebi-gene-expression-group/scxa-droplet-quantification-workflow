#!/usr/bin/env bash

usage() { echo "Usage: $0 [-r <gzipped read file>] [-u UMI length] [-b cell barcode length] [-n number of reads to test]" 1>&2; }

r=
u=
b=
n=1000

while getopts ":r:u:b:n:" o; do
    case "${o}" in
        r)
            r=${OPTARG}
            ;;
        u)
            u=${OPTARG}
            ;;
        b)
            b=${OPTARG}
            ;;
        n)
            n=${OPTARG}
            ;;
        *)
            usage
            exit 0
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${r}" ] || [ -z "${u}" ] || [ -z "${b}" ] || [ -z "${n}" ]; then
    usage
    exit 1
fi

if [ ! -e "$r" ]; then
    echo "Read file $r does not exist" 1>&2
    exit 1
fi

# Calculate lenths of first n reads

echo "Checking $n reads"

lengths=$(zcat $r | \
    head -n $((4 * $n)) | \
    sed -n '2~4p' | \
    awk '{print length()}' | \
    sort -nr | uniq)

# Find number of unque lengths

nLengths=$(echo -e "$lengths" | wc -l)
if [ "$nLengths" == 1 ]; then
    qualifier='all'
else
    qualifier='some'
fi

# Find longest read

longest=$(echo -e "$lengths" | head -n 1)

# Reads should be at least the UMI + CB length

targetLength=$(($u + $b))

# Print warnings for fishy things, die for big issues

if [ "$nLengths" -gt 1 ]; then
    echo "WARNING: UMI/ barcode reads are of variable length" 1>&2
fi

if [ "$longest" -gt "$targetLength" ]; then
    echo "WARNING: $qualifier UMI/ barcode reads in $(basename $r) are longer than UMI + cell barcode (=$targetLength), max length: $longest. This happens with sequencer run-on, but make sure you're sure of the barcode configuration" 1>&2
elif [ "$longest" -lt "$targetLength" ]; then
    echo "[ERROR} No reads in $(basename $r) meet UMI + cell barcode length threshold of $targetLength (max length: $longest), these reads cannot pass to droplet quantification." 1>&2
    exit 1
else
    echo "[SUCCESS: $qualifier reads in $(basename $r) match UMI/ cell barcode length threshold of $targetLength (max length: $longest)"
fi
