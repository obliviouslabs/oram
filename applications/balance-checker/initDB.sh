export LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH
DBPath=db_usdt
BalancePath=usdt_balance.txt
rm -r $DBPath
g++ init_db_balance.cpp -o init_db_balance -L/usr/local/lib -lrocksdb
./init_db_balance $DBPath $BalancePath