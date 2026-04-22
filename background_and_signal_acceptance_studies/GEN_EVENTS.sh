for model in momentum_constrained core floating; do

    #python generate_signal_events_CL.py --masses-a 0.22 --masses-dm 1000 10000 100000 500000 \
                    python generate_signal_events_CL.py --masses-a 0.22 --masses-dm 10000 20000 30000 40000 50000 60000 70000 80000 90000 \
                                        --depth -8 2000 --disk-radius 2000 \
                                        --save-to-file \
                                        --average-eloss \
                                        --output-directory OUTPUT_FILES \
                                        --dm-model $model  \
                                        --ntrials 2000 \
                                        --nevents 1_000_000

done

#python generate_signal_events_CL.py --masses-a 0.22 --masses-dm 1000 10000 100000 500000 \
                                    #--depth -10 \
                                    #--save-to-file \
                                    #--output-directory OUTPUT_FILES \
                                    #--average-eloss \
                                    #--dm-model momentum_constrained  \
                                    #--ntrials 2000 \
                                    #--nevents 1_000_000
                                    #--ntrials 20 \
                                    #--nevents 1_000_000 

#python generate_signal_events_CL.py --masses-a 0.22 --masses-dm 20000 \
                                    #--depth -8 \
                                    #--save-to-file \
                                    #--output-directory OUTPUT_FILES \
                                    #--average-eloss \
                                    #--dm-model momentum_constrained  \
                                    #--ntrials 2 \
                                    #--nevents 1_000_000
                                    #--ntrials 20 \
                                    #--nevents 1_000_000 

#python generate_signal_events_CL.py --masses-a 0.22 --masses-dm 20000 \
                                    #--nevents 1_000_000 --ntrials 20 \
                                    #--depth -8 -6000 --disk-radius 40 \
                                    #--dm-model core  \
                                    #--save-to-file \
                                    #--output-directory OUTPUT_FILES
                                    #--average-eloss \
                                    #--dm-model momentum_constrained  \

#python generate_signal_events_CL.py --masses-a 0.22 --masses-dm 7000 \
                                    #--nevents 1_000_000 --ntrials 5000 \
                                    #--depth -8 -3500 --disk-radius 3500 \
                                    #--dm-model floating  \
                                    #--save-to-file \
                                    #--output-directory OUTPUT_FILES
#
#python generate_signal_events_CL.py --masses-a 0.22 --masses-dm 7000 \
                                    #--nevents 1_000_000 --ntrials 5000 \
                                    #--depth -8 -3500 --disk-radius 3500 \
                                    #--dm-model momentum_constrained  \
                                    #--save-to-file \
                                    #--output-directory OUTPUT_FILES
#
#python generate_signal_events_CL.py --masses-a 0.22 --masses-dm 7000 \
                                    #--nevents 1_000_000 --ntrials 5000 \
                                    #--depth -8 -3500 --disk-radius 3500 \
                                    #--dm-model floating  \
                                    #--save-to-file \
                                    #--output-directory OUTPUT_FILES \
                                    #--average-eloss 
#
#python generate_signal_events_CL.py --masses-a 0.22 --masses-dm 7000 \
                                    #--nevents 1_000_000 --ntrials 5000 \
                                    #--depth -8 -3500 --disk-radius 3500 \
                                    #--dm-model momentum_constrained  \
                                    #--save-to-file \
                                    #--output-directory OUTPUT_FILES \
                                    #--average-eloss 
