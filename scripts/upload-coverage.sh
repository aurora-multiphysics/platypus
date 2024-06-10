#!/bin/bash

cd /opt/platypus || exit 1

for filename in build/**/*.C; do
    gcov -n -o . "$filename" > /dev/null
done

bash <(curl -s https://codecov.io/bash)
