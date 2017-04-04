#!/bin/bash

if [[ -v $CONDA_PY ]]; then
  # conda build
  if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
      brew install md5sha1sum
  fi
  source devtools/travis-ci/install_miniconda.sh
  conda install --yes --file requirements.txt
else
  # plain build
  pip install -r requirements.txt
fi
