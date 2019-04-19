#!/bin/bash
# rm old archives
rm dist/ build/ isoCirc.egg-info/ -rf 2> /dev/null
# generate distribution archives
python setup.py sdist bdist_wheel
# upload package
twine upload dist/*

