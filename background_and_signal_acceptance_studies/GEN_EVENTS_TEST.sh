NEVENTS=1_000
NTRIALS=5

for model in "floating"; do

    echo $model
    echo $NEVENTS
    echo $NTRIALS

    python generate_signal_events_CL_v2.py \
        --masses-a 0.22 --masses-dm 10 20 30 40 50 60 70 80 90  \
        --depth -8 -100 \
        --disk-radius 100 \
        --save-to-file --average-eloss \
        --output-directory OUTPUT_FILES \
        --data-root data \
        --dm-model $model \
        --ntrials $NTRIALS --nevents $NEVENTS


done
