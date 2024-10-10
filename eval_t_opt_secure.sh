max=5
for i in `seq 0 $max`
do
    for j in `seq $((i+1)) $max`
    do
        echo "Genomes" $i $j
        python3 edit_dist_runner.py $i $j 2041
        python3 dataProcessing/generate_player_data.py 0 $i $j
        cp Player-Data/Input-P0-0 ../MP-SPDZ/Player-Data/Input-P0-0
        cp Player-Data/Input-P1-0 ../MP-SPDZ/Player-Data/Input-P1-0
        cd ..
        cd MP-SPDZ
        ./yao-party.x -p 0 runner_phase1_2k_48
        cd ..
        cd secure-edit-distance
    done
done