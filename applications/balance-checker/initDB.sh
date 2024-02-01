DBPath=db_rcc
BalancePath=rcc_balance.txt
rm -r $DBPath
g++ init_db_balance.cpp -o init_db_balance -L/usr/local/lib -lrocksdb
./init_db_balance $DBPath $BalancePath