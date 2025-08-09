#for energy in 50 100 200 300 400 500 600 700 800 900 1000; do
#for energy in 2000 3000 4000 5000 6000 7000 8000 9000 10000; do
for energy in 20000 30000 40000 50000 60000 70000 80000 90000 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000; do
#for energy in 50 100; do
    OUTMAC=run_"$energy".mac
    OUTCSV=muonSteps_"$energy"_GeV.csv
    cat run_EMPTY.mac  | sed s/XXXENERGYXXX/$energy/ > $OUTMAC
    echo ./build/MuonThroughRock $OUTMAC
         ./build/MuonThroughRock $OUTMAC
    echo mv muonSteps.csv $OUTCSV
         mv muonSteps.csv $OUTCSV
    echo
done
