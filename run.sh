cd src
echo "Recompiling any new changes made to program!"
make
cd ..
echo "Estimation has begun!"
./CyGMM > out.txt
echo "Estimation finished! Please look in out.txt for newly computed estimation results!"