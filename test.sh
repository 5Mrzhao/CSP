input="/home/zyx/Circuitmatrix/mtx_set.csv"
{
  read
  i=1
  while IFS=',' read -r Name
  do
    echo "$Name"
    ./main $Name
    i=`expr $i + 1`
  done 
} < "$input"