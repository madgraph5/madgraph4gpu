date
while IFS= read -r line
do
  cmd="./check.exe -p $line"
  echo
  echo $cmd
  $cmd
  date
done < $1
