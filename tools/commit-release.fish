#!/usr/bin/env fish

if type -q yq
    set VERSION (yq .version CITATION.cff)

    git add LICENSE.txt
    git add CITATION.cff
    git add Project.toml

    git commit -m "Change version to $VERSION"

    echo "Version $VERSION committed (but not pushed)"
else
    echo "Install yq (from https://github.com/mikefarah/yq)"
end
