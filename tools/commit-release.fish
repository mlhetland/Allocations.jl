#!/usr/local/bin/fish

set VERSION (yq .version CITATION.cff)

git add LICENSE.txt
git add CITATION.cff

git commit -m "Change version to $VERSION"
git tag -a "v$VERSION" -m "Version $VERSION"

echo "Version $VERSION committed and tagged (but not pushed)"
