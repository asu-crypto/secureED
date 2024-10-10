cd ..
cd MP-SPDZ
max=5
for i in `seq 0 $max`
do
    for j in `seq $((i+1)) $max`
    do
        ./yao-party.x -p 1 runner_phase1_2k_48
    done
done