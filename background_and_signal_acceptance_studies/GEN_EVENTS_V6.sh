NEVENTS=1_000_000
NTRIALS=2000

for model in "floating"; do

    echo $model
    echo $NEVENTS
    echo $NTRIALS

    python generate_signal_events_CL_v2.py \
        --masses-a 0.22 --masses-dm 10000 20000 30000 40000 50000 60000 70000 80000 90000  \
        --depth -8 -6000 \
        --disk-radius 6000 \
        --save-to-file --average-eloss \
        --output-directory OUTPUT_FILES \
        --data-root data \
        --dm-model $model \
        --ntrials $NTRIALS --nevents $NEVENTS
        #--merge-run-id 20260611T02503449b0
        #--seed 1000 \
        #--dm-model momentum_constrained \
        #--masses-a 0.22 --masses-dm 10 20 30 40 50 60 70 80 90  \


done
