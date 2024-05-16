#!/usr/bin/env fish

if type -q yq
    set VERSION (yq .version CITATION.cff)
    git tag -a "v$VERSION" -m "Version $VERSION"

    echo "Version $VERSION tagged (but not pushed)"
else
    echo "Install yq (from https://github.com/mikefarah/yq)"
end
