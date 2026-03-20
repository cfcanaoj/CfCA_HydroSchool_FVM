# Brio-Wu衝撃波管問題

このディレクトリには，[Brio & Wu (1988)](https://ui.adsabs.harvard.edu/link_gateway/1988JCoPh..75..400B/doi:10.1016/0021-9991(88)90120-9)が提案した衝撃波管問題を解くサンプルコードと結果を可視化するためのファイルが入っている。

## 初期条件

初期不連続面の左状態と右状態は以下のように設定する。

$$
\begin{pmatrix}
\rho_{\rm L} \\
v_{x,{\rm L}} \\
v_{y,{\rm L}} \\
v_{z,{\rm L}} \\
P_{\rm L} \\
B_{x,{\rm L}} \\
B_{y,{\rm L}} \\
B_{z,{\rm L}} 
\end{pmatrix}
=\begin{pmatrix}
1 \\
0 \\
0 \\
0 \\
1 \\
0.75 \\
1 \\
0 \\
\end{pmatrix},
\qquad
\begin{pmatrix}
\rho_{\rm R} \\
v_{x,{\rm R}} \\
v_{y,{\rm R}} \\
v_{z,{\rm R}} \\
P_{\rm R} \\
B_{x,{\rm R}} \\
B_{y,{\rm R}} \\
B_{z,{\rm R}} 
\end{pmatrix}
=\begin{pmatrix}
0.125 \\
0 \\
0 \\
0 \\
0.1 \\
0.75 \\
-1 \\
0 \\
\end{pmatrix},
$$

複合波を含む厳密解は以下の図のようになる。
<p align="center">
  <img src="./briowu_nonregsol.png" width="600">
</p>

## ディレクトリ内の構造

ファイル名など     | 説明  
------------------|----------
 main.f90         | 流体学校の実習用。複数のサブルーチンが未完
 hlld.f90         | HLLD流束を求めるサブルーチン
 ans/main.f90     | 完成品
 briowu_nonregsol.dat | Brio-Wu衝撃波管問題の複合波を含む厳密解
 briowu_nonregsol.png | Brio-Wu衝撃波管問題の複合波を含む厳密解を示した図
 可視化スクリプト   | RealtimeAnim.plt, MakeAnime.sh, MakeAnime.py, MakeCompare.py


## 可視化

### 動画作成

**RealtimeAnim.plt**

gnuplotを使ってプロットを連続的に表示し，動画として表示するためのスクリプト。
dirname内のibeg番からifin番までの動画を表示するには，以下を実行する。

```bash
gnuplot -c RealtimeAnim.plt ibeg ifin dirname
```

引数を付けない場合
```bash
gnuplot RealtimeAnim.plt
```
は，ibeg=0, ifin=40, dirname=laxが採用される。

**MakeAnime.sh**

gnuplotを使ってdirname内のibeg番からifin番までの画像ファイルを作成し，
ffmpegを使って動画ファイルを作成するためのスクリプト。

```bash
./MakeAnime.sh ibeg ifin dirname
```

**MakeAnim.py**

pythonを使って，特定の範囲のスナップショットから動画を作成するスクリプト。

Usage: python3 MakeAnime.py [step_s] [step_e] [dir1] [dir2] ...

```bash
python3 MakeAnime.py
```

laxとhdc内のsnap00000.datからsnap00040.datまでのデータから動画を作成する場合は，以下を実行する。
```bash
python3 MakeAnime.py 0 40 lax hdc
```

### スナップショットの図示

**MakeCompare.py**

pythonを使って，複数の出力ディレクトリ内の特定snapXXXXX.datを比較するためのスクリプト

Usage: python MakeAnime.py [step_s] [step_e] [dir1] [dir2] ... 

lax/snap00040.datとhll/snap00040.datを比較する場合は，以下を実行する。
```bash
python3 MakeAnime.py 40 lax hdc
```
