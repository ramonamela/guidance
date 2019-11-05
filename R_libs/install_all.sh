#!/bin/bash

install_dir=${1}

R CMD INSTALL *tar.gz --library=${install_dir}

