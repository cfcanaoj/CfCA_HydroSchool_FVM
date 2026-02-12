#!/bin/bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 dir nstart nend" >&2
  exit 1
fi

dir=$1
nmin=$2
nmax=$3
prename=E

if [ ! -d "$dir" ]; then
  echo "cannot find $dir" >&2
  exit 1
fi

for n in $(seq "$nmin" "$nmax"); do
  # dir は gnuplot に「文字列」として渡す（シングルクォートで囲む）
  gnuplot -c  MakePngFile.plt "${dir}" "${n}"
done

fstfile=$(ls -1 "${dir}/${prename}"*.png 2>/dev/null | head -1 || true)
if [ -z "$fstfile" ]; then
  echo "no png files found in $dir" >&2
  exit 1
fi

fstnum=$(basename "$fstfile" | tr -cd '0-9\n' | sed -e 's/^0\+\([0-9]\+\)$/\1/')
echo "$fstnum"

echo "wmv format"
ffmpeg -y -r 10 -start_number "$fstnum" -i "${dir}/${prename}%5d.png" -b 6000k -vcodec wmv2 -pass 1 -r 10 -an "${dir}/animation.wmv"

echo "mp4 format"
ffmpeg -y -r 10 -start_number "$fstnum" -i "${dir}/${prename}%5d.png" -vcodec libx264 -pix_fmt yuv420p -r 10 -an "${dir}/animation.mp4"
