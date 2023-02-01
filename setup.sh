#!/bin/bash
SETUP_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export IGMF_ROOT=$SETUP_DIR
echo "IGMF_ROOT set to " $IGMF_ROOT
export PYTHONPATH=${PYTHONPATH}:$IGMF_ROOT
echo "PYTHONPATH set to " $PYTHONPATH
export PATH=$IGMF_ROOT/bin:$PATH
echo "PATH set to " $PATH
