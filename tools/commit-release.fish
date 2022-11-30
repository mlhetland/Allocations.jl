#!/usr/bin/env fish

if type -q yq
    set VERSION (yq .version CITATION.cff)

    git add LICENSE.txt
    git add CITATION.cff
    git add Project.toml

    git commit -m "Change version to $VERSION"
    git tag -a "v$VERSION" -m "Version $VERSION"

    echo "Version $VERSION committed and tagged (but not pushed)"
else
    echo "Install yq (from https://github.com/mikefarah/yq)"
end
