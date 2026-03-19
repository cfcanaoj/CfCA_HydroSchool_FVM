# 音波の伝播

このディレクトリには，音波の伝播を解くサンプルコードと結果を可視化するためのファイルが入っている。

## 初期条件

非摂動状態として，静止した一様ガスを考える($\rho=1, v_x = 0, P=1/\gamma$)。

$$ 
\rho = A \sin (kx)
$$

$\gamma=1.4$を採用している。

## コンパイル

```bash
gfortran -O3 main.f90
```

## 可視化

* RealtimeAnim.plt

gnuplotを使って画像を連続的に表示し，動画として表示するためのスクリプト。
dirname内のibeg番からifin番までの動画を表示するには，以下を実行する。

```bash
gnuplot -c RealtimeAnim.plt ibeg ifin dirname
```

* MakeAnime.sh

gnuplotを使ってdirname内のibeg番からifin番までの画像ファイルを作成し，
ffmpegを使って動画ファイルを作成するためのスクリプト。

```bash
./MakeAnime.sh ibeg ifin dirname
```


* MakeCompare.py

pythonを使って，複数の出力ディレクトリ内のスナップショットを比較するためのプログラム

Usage: python MakeAnime.py [step_s] [step_e] [dir1] [dir2] ... 

* MakeAnim.py

pythonを使って，複数の出力ディレクトリ内のスナップショットを比較するためのプログラム

Usage: python3 MakeCompare.py [step] [dir1] [dir2] ...
