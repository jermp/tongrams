#!/usr/bin/env sh

python3 -m build -nwx .
sudo python3 -m pip install --upgrade ./dist/*.whl
