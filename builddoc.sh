#!/bin/bash

set -eu

# Build docs
julia docs/make.jl

# Switch Branch
git checkout gh-pages

# Copy files
git checkout master docs/build/


# Relocate files
cp -r docs/build/. docs/

# commit all changes

# git add docs/
# git commit -m "building docs"

# git push origin gh-pages

# Switch Branch
git checkout master
