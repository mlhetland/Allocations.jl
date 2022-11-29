#!/usr/local/bin/fish

set YEAR (date +"%Y")
set DATE (date +"%Y-%m-%d")

set NAMES (yq '.authors | map(.given-names + " " + .family-names)
              | join(", ")' CITATION.cff)

read -P "Version: " VERSION

# https://github.com/mikefarah/yq
yq -i ".version = \"$VERSION\" | .date-released = \"$DATE\"" CITATION.cff

set CPREFIX "Copyright (c)"
set COPYRIGHT (echo "$CPREFIX 2020-$YEAR $NAMES" | fold -w 78)

sed -e "s/^$CPREFIX.*/$COPYRIGHT/" LICENSE.txt > LICENSE.updated
if cmp -s LICENSE.txt LICENSE.updated
    rm LICENSE.updated
    echo "CITATION.cff updated"
else
    mv LICENSE.updated LICENSE.txt
    echo "LICENSE.txt and CITATION.cff updated"
end

