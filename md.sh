#!/bin/bash


pid=-1

run () {
    cat $2 | ./$1
    echo $!;
    return $!;
}

kill () {
    kill $1;
}

usage () {
    echo "Usage: ./md.sh <MD> <inputfile>"
}

attach () {
    pid=$1;
    tail -f /proc/$pid/fd/1;
}

if [ $# -ne 2 ]; then
    usage;
    exit 1;
fi

pid="$(run $1 $2)";
echo $pid;