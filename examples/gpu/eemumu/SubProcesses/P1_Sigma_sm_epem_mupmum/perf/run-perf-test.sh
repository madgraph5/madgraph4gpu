date
while IFS= read -r line
do
  sleep 30
  cmd="./check.exe -p $line"
  echo
  echo $cmd
  $cmd
  date
done < $1
