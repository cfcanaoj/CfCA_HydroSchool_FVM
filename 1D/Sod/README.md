# Sod Shock Tube Test

このディレクトリには，衝撃波管問題 Sod解を解くサンプルコードと結果を可視化するためのファイルが入っている。

## 初期条件

初期不連続面の左状態と右状態は以下のように設定する。

$$ {\left(\matrix{\rho_{\rm L} \cr v_{\rm L} \cr P_{\rm L} } \right)}
= {\left(\matrix{1 \cr 0 \cr 1 } \right)},
{\left(\matrix{\rho_{\rm R} \cr v_{\rm R} \cr P_{\rm R} } \right)}
= {\left(\matrix{0.125 \cr 0 \cr 0.1 } \right)}
$$

$\gamma=1.4$を採用している。

## コンパイル

gfortran -O3 main.f90

## 可視化

* RealtimeAnim.plt

gnuplotを使って画像を連続的に表示し，動画として表示するためのスクリプト。
dirname内のibeg番からifin番までの動画を表示するには，以下を実行する。

`gnuplot -c RealtimeAnim.plt ibeg ifin dirname`

* MakeAnime.sh

gnuplotを使ってdirname内のibeg番からifin番までの画像ファイルを作成し，
ffmpegを使って動画ファイルを作成するためのスクリプト。

`./MakeAnime.sh ibeg ifin dirname`


* MakeCompare.py

pythonを使って，複数の出力ディレクトリ内のスナップショットを比較するためのプログラム

Usage: python MakeAnime.py [step_s] [step_e] [dir1] [dir2] ... 

* MakeAnim.py

pythonを使って，複数の出力ディレクトリ内のスナップショットを比較するためのプログラム

Usage: python3 MakeCompare.py [step] [dir1] [dir2] ...
