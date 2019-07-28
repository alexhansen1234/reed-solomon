err=0;
tests=0;
for i in {1..100};
do
  ../reed-solomon 1>/dev/null;
  if [ $? != 0 ];
  then
      echo Test $tests failed.
      err=$(($err + 1));
  else
      echo Test $tests passed.
  fi
  tests=$(($tests+1));
done;
echo "Total errors $err in $tests tests."
