#!/bin/bash

case "$1" in
    gc)
        python3 dataProcessing/generate_player_data.py 0 "$2" "$3";;
    ss)
        python3 dataProcessing/generate_player_data.py 1 "$2" "$3";;
esac
    cp Player-Data/Input-P0-0 ../MP-SPDZ/Player-Data/Input-P0-0
    cp Player-Data/Input-P1-0 ../MP-SPDZ/Player-Data/Input-P1-0