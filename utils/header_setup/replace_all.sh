#!/bin/bash

  SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

  # Change to the GUIDANCE root directory
  cd "$SCRIPT_DIR"/../../ || exit 1

  # Add java headers
  find . -name "*.java" -exec "$SCRIPT_DIR"/replace_header.sh {} \; 

