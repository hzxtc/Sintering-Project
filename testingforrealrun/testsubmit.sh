#! /bin/bash
# BP case submit file

now1=$(date +"%T")
echo "Begin time : $now1"
echo "1: "
python3 sintering_BP.py
tar -czvf metropolis1.tar.gz metropolis
echo "2: "
python3 sintering_BP.py
tar -czvf metropolis2.tar.gz metropolis
echo "3: "
python3 sintering_BP.py
tar -czvf metropolis3.tar.gz metropolis
echo "4: "
python3 sintering_BP.py
tar -czvf metropolis4.tar.gz metropolis
echo "5: "
python3 sintering_BP.py
tar -czvf metropolis5.tar.gz metropolis

now2=$(date +"%T")
echo "FINISHED time : $now2"
