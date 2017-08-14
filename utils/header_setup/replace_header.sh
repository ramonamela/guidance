#!/bin/bash

  SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

  # Check and get parameters
  if [ $# != 1 ]; then
    echo "Usage: replace_header.sh <file>"
    exit 1
  fi

  file=$1

  # Replace header
  cat "$SCRIPT_DIR"/header.template.java > "$SCRIPT_DIR"/tmp
  awk -f "$SCRIPT_DIR"/remove_header_java.awk < "$file" >> "$SCRIPT_DIR"/tmp
  mv "$SCRIPT_DIR"/tmp "$file"

