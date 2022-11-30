#!/usr/bin/env fish

# https://github.com/mikefarah/yq
set VERSION (yq .version CITATION.cff)

git add LICENSE.txt
git add CITATION.cff
git add Project.toml

git commit -m "Change version to $VERSION"
git tag -a "v$VERSION" -m "Version $VERSION"

echo "Version $VERSION committed and tagged (but not pushed)"
