#!/bin/bash

#declare -i nmin nmax

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 dir nstart nend ."
  echo "Error: at least 3 arguments are required." >&2
  exit 1
fi

dir=$1
nmin=$2
nmax=$3
prename=E

if [ ! -d $dir ]; then
    echo "cannot find "$dir
    exit
fi

for n in `seq $nmin $nmax` ;do
gnuplot -e "num=$n; dir=$dir" MakePngFile.plt
done

# first file
fstfile=`ls -1 ${dir}/${prename}*.png  2>/dev/null | head -1`
declare -i fstnum=`echo  ${fstfile##*/} | tr -cd '0123456789\n' |sed -e 's/^0\+\([0-9]\+\)$/\1/'`
echo $fstnum

echo "wmv format"
ffmpeg -y -r 10  -start_number ${fstnum} -i ${dir}/${prename}%5d.png -b 6000k -vcodec wmv2 -pass 1 -r 10 -an ${dir}/animate.wmv

echo "mp4 format"
ffmpeg -y -r 10  -start_number ${fstnum} -i ${dir}/${prename}%5d.png -vcodec libx264 -pix_fmt yuv420p -r 10 -an ${dir}/animate.mp4


exit
