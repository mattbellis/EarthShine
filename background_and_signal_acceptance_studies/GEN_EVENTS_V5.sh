NEVENTS=1_000_000
NTRIALS=500

for model in floating; do

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
        #--merge-run-id 20260611T02503449b0
        #--seed 1000 \
        #--dm-model momentum_constrained \
        #--masses-a 0.22 --masses-dm 10 20 30 40 50 60 70 80 90  \


    python generate_signal_events_CL_v2.py \
        --masses-a 0.22 --masses-dm 100 200 300 400 500 600 700 800 900  \
        --depth -8 -1000 \
        --disk-radius 1000 \
        --save-to-file --average-eloss \
        --output-directory OUTPUT_FILES \
        --data-root data \
        --dm-model $model \
        --ntrials $NTRIALS --nevents $NEVENTS
        #--merge-run-id 20260611T02503449b0
        #--seed 1000 \
        #--dm-model momentum_constrained \
        #--masses-a 0.22 --masses-dm 10 20 30 40 50 60 70 80 90  \

    python generate_signal_events_CL_v2.py \
        --masses-a 0.22 --masses-dm 1000 2000 3000 4000 5000 6000 7000 8000 9000  \
        --depth -8 -4000 \
        --disk-radius 4000 \
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

NEVENTS=1_000_000
NTRIALS=10

for model in core; do

    echo $model
    echo $NEVENTS
    echo $NTRIALS

    python generate_signal_events_CL_v2.py \
        --masses-a 0.22 --masses-dm 10 20 30 40 50 60 70 80 90  \
        --depth -8 -100 \
        --disk-radius 40 \
        --save-to-file --average-eloss \
        --output-directory OUTPUT_FILES \
        --data-root data \
        --dm-model $model \
        --ntrials $NTRIALS --nevents $NEVENTS
        #--merge-run-id 20260611T02503449b0
        #--seed 1000 \
        #--dm-model momentum_constrained \
        #--masses-a 0.22 --masses-dm 10 20 30 40 50 60 70 80 90  \


    python generate_signal_events_CL_v2.py \
        --masses-a 0.22 --masses-dm 100 200 300 400 500 600 700 800 900  \
        --depth -8 -1000 \
        --disk-radius 40 \
        --save-to-file --average-eloss \
        --output-directory OUTPUT_FILES \
        --data-root data \
        --dm-model $model \
        --ntrials $NTRIALS --nevents $NEVENTS
        #--merge-run-id 20260611T02503449b0
        #--seed 1000 \
        #--dm-model momentum_constrained \
        #--masses-a 0.22 --masses-dm 10 20 30 40 50 60 70 80 90  \

    python generate_signal_events_CL_v2.py \
        --masses-a 0.22 --masses-dm 1000 2000 3000 4000 5000 6000 7000 8000 9000  \
        --depth -8 -4000 \
        --disk-radius 40 \
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
