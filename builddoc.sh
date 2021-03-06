#!/bin/bash

set -eu

# Build docs
cd docs
make
cd ..

# Copy files
mkdir ../temp
cp -r docs/build/. ../temp/

# Switch Branch
git checkout gh-pages

# Relocate files
cp -r ../temp/. docs/

# commit all changes

git add docs/
git commit -m "building docs"

# git push origin gh-pages

# Switch Branch
git checkout master

rm -rf ../temp/
