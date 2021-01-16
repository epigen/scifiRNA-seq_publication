


# I copied the CSV files from my hard drive backup I made when I left Vienna
# into a folder called 'csvs'

# This script will reconstruct the original structure needed to run the scripts
# by making sample-wise directories and softlinking the files to them

for F in `ls csvs/*metrics.csv.gz`; do
    echo $F

    DIR=`basename ${F/.metrics.csv.gz/}`
    mkdir -p data/${DIR}
    cd data/${DIR}
    ln -s ../../$F ./
    cd -
done
