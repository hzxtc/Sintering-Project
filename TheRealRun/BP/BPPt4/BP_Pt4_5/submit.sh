#! /bin/bash
# BP case submit file

echo "-----------------------------------------------------------------"
now1=$(date +"%T")
echo "Begin time : $now1"
echo "1: "
python3 sintering_BP.py
tar -czvf metropolis1.tar.gz metropolis
now2=$(date +"%T")
echo "time : $now2"
echo "2: "
python3 sintering_BP.py
tar -czvf metropolis2.tar.gz metropolis
now3=$(date +"%T")
echo "time : $now3"
echo "3: "
python3 sintering_BP.py
tar -czvf metropolis3.tar.gz metropolis
now4=$(date +"%T")
echo "time : $now4"
echo "4: "
python3 sintering_BP.py
tar -czvf metropolis4.tar.gz metropolis
now5=$(date +"%T")
echo "time : $now5"
echo "5: "
python3 sintering_BP.py
tar -czvf metropolis5.tar.gz metropolis

now6=$(date +"%T")
echo "FINISHED time : $now6"
echo "-----------------------------------------------------------------"
