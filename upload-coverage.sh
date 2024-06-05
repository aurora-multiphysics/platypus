#!/bin/bash

cd /opt/platypus

for filename in `find src -name "*.C"`; 
do 
  gcov -n -o . $filename > /dev/null; 
done

bash <(curl -s https://codecov.io/bash)
