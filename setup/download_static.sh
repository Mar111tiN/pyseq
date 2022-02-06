#!/bin/sh

# downloads the STATIC data for use with genomic tools

STATICTXT="https://www.dropbox.com/s/urmw86jw3xv7mtb/static_links.txt"

# store current folder
CURRENT=$(pwd);
if [ -z "$1" ]; then
    echo "Please provide a static folder as argument!"
else
    FOLDER="$1";
    echo "Downloading static data to $FOLDER";
    cd $FOLDER;
    # download the links txt
    wget $STATICTXT;
    # download the first files from the link list
    cat static_links.txt | awk 'NR<3' | xargs wget;

    # check and expand static folder
    STATIC="hg38_static.tar.gz"
    echo "Expanding static tar.gz files and checking integrity";
    md5sum -c ${STATIC}.md5 && tar -xzvf $STATIC && rm ${STATIC}* && echo "Static data downloaded to $FOLDER";
    # rollback
    cd $FOLDER && rm static_links.txt;
    cd $CURRENT;
fi