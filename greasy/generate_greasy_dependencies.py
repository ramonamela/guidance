#!/usr/bin/env python

import sys




if __name__ == "__main__":
    input_file = sys.argv[1]
    with open(input_file) as infile:
            lines = (line.rstrip() for line in infile)  # All lines including the blank ones
            lines = (line for line in lines if not line[0] == "#" and line.strip())  # Non-blank lines
            lines = (map(lambda x: x.strip(), line.split("=")) for line in lines)
            param_dict = dict(lines)
            for line in lines:
                print(line + " " + lines[line])
