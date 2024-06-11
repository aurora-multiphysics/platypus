#!/bin/bash

cd /opt/platypus || exit 1
bash <(curl -s https://codecov.io/bash)
