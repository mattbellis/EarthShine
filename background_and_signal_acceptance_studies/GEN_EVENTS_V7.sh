for model in "floating"; do

    NEVENTS=1_000_000
    NTRIALS=1000

    echo $model
    echo $NEVENTS
    echo $NTRIALS

    python generate_signal_events_CL_v2.py \
        --masses-a 0.22 --masses-dm 100 \
        --depth -8 -4000 \
        --disk-radius 4000 \
        --save-to-file --average-eloss \
        --output-directory OUTPUT_FILES \
        --data-root data \
        --dm-model $model \
        --ntrials $NTRIALS --nevents $NEVENTS

done

for model in "floating"; do

    NEVENTS=1_000_000
    NTRIALS=1000

    echo $model
    echo $NEVENTS
    echo $NTRIALS

    python generate_signal_events_CL_v2.py \
        --masses-a 0.22 --masses-dm 100 \
        --depth -8 -4000 \
        --disk-radius 4000 \
        --save-to-file --average-eloss \
        --output-directory OUTPUT_FILES \
        --data-root data \
        --dm-model $model \
        --ntrials $NTRIALS --nevents $NEVENTS

done

for model in "floating"; do

    NEVENTS=1_000_000
    NTRIALS=1000

    echo $model
    echo $NEVENTS
    echo $NTRIALS

    python generate_signal_events_CL_v2.py \
        --masses-a 0.22 --masses-dm 100 \
        --depth -8 -4000 \
        --disk-radius 4000 \
        --save-to-file --average-eloss \
        --output-directory OUTPUT_FILES \
        --data-root data \
        --dm-model $model \
        --ntrials $NTRIALS --nevents $NEVENTS

done
