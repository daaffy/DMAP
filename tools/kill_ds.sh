#!/bin/bash

# https://stackoverflow.com/questions/107701/how-can-i-remove-ds-store-files-from-a-git-repository
find . -name .DS_Store -print0 | xargs -0 git rm -f --ignore-unmatch