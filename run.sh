cd src
echo "Recompiling any new changes made to program!"
make
echo "Estimation has begun!"
./sig > ../out.txt
echo "Estimation finished! Please look in out.txt for newly computed estimation results!"
cd ..