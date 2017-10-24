#!/bin/bash

#PRESET="veryslow"
#PRESET="slow"
PRESET="medium"
#PRESET="ultrafast"

CRF=21

TYPE="mp4"

if [ "$4" = "" ]; then
    echo "USAGE: `basename $0` <begining_files> <format (jpg,png,...)> <height video (pixels)> <fps>"
    exit
fi

SCALE="-vf scale='-2:${3}'"
FPS=${4}

#Codec stuff

if [ "${TYPE}" = "mp4" ]; then
    VC="-c:v libx264 -profile:v high444"
    info="x264_${3}px"
elif [ "${TYPE}" = "webm" ]; then
    VC="-c:v libvpx"
    info="vpx_${3}px"
elif [ "${TYPE}" = "mov" ]; then
    VC="-c:v h264_nvenc"
    info="h264_${3}px"
else
    echo "Boo!" ; exit
fi

fo="movie_${1}_${info}.${TYPE}"

rm -f ${fo}

ffmpeg -f image2 -framerate ${FPS} \
       -pattern_type glob -i "${1}*.${2}" \
       ${VC} -preset ${PRESET} \
       -crf ${CRF} -refs 16 ${SCALE} \
       -pix_fmt yuv420p \
       ${fo}

echo ; echo ; echo
exit 0
