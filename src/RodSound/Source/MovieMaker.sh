#!/bin/sh

#  MovieMaker.sh
#  Visualizer
#
#  Created by eschweickart on 6/12/14.
#

# Move to the result directory
if [ ! -d $1 ]; then
  echo "Error: Result directory does not exist";
  exit 1;
fi
pushd $1;

# Encode video using ffmpeg
# 60 fps, mpeg2 encoding
# Overwrites existing files of the same name
/usr/local/bin/ffmpeg -y -r 60 -i png/%d.png -i ./result.wav -c:v mpeg2video -q:v 5 ./$2.mpg;

# TODO: delete pngs?


popd;