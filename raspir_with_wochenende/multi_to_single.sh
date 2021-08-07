#!/bin/bash

# Marie-Madlen Pust
# Merge all files and replace header with file name

awk -F "," '
    BEGIN { header = "organism" }
    {
        key = $1
        values[key] = values[key] FS $2
    }
    FNR == 1 {
        a = split(FILENAME,array, "_")
        header = header FS array[1]"_"array[2]
    }
    END {
        print header
        for (key in values)
            print key values[key]
    }
' *.wochenende.rep.us.raspir.csv > wochenende.rep.us.raspir.merged.csv

sed -i '/  ,,,/d' wochenende.rep.us.raspir.merged.csv
