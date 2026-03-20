# Sod衝撃波管問題

このディレクトリには，衝撃波管問題 Sod解を解くサンプルコードと結果を可視化するためのファイルが入っている。

## 初期条件

初期不連続面の左状態と右状態は以下のように設定する。

$$
\begin{pmatrix}
\rho_{\rm L} \\
v_{\rm L} \\
P_{\rm L}
\end{pmatrix}
=\begin{pmatrix}
1 \\
0 \\
1
\end{pmatrix},
\qquad
\begin{pmatrix}
\rho_{\rm R} \\
v_{\rm R} \\
P_{\rm R}
\end{pmatrix}=
\begin{pmatrix}
0.125 \\
0 \\
0.1
\end{pmatrix}
$$

$\gamma=1.4$を採用している。

## ディレクトリ内の構造

ファイル名         | 説明  
------------------|----------
 main.f90         | 流体学校の実習用。複数のサブルーチンが未完
 hllc.f90         | hllc流束を求めるサブルーチン
 ans/main.f90     | 完成品
 RiemannSolver.py | 任意の初期条件からRiemann問題の解を求めるスクリプト
 可視化スクリプト   | RealtimeAnim.plt, MakeAnime.sh, MakeAnime.py, MakeCompare.py


## Riemann問題を求解

以下を実行すると，Sod解の厳密解がファイル(riemann_exact.dat)として出力される。
```bash
>python3 RiemannSolver.py 

=== adopted initial states ===
rhoL  = 1.000000   | rhoR  = 0.125000
velL  = 0.000000   | velR  = 0.000000
preL  = 1.000000   | preR  = 0.100000

gamma = 1.400000
time  = 0.200000
x0    = 0.000000
outputfile = riemann_exact.dat
```
ファイル名を変更するには，オプション「-o」で指定する。たとえば，sod_ana.datとする場合，
```bash
>python3 RiemannSolver.py -o sod_ana.dat
```
を実行する。

他のオプションについては，
```bash
>python3 RiemannSolver.py -h
```
を参照のこと。



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

* MakeCompare.py

pythonを使って，複数の出力ディレクトリ内の特定snapXXXXX.datを比較するためのスクリプト

Usage: python MakeAnime.py [step_s] [step_e] [dir1] [dir2] ... 

lax/snap00040.datとhll/snap00040.datを比較する場合は，以下を実行する。
```bash
python3 MakeAnime.py 40 lax hdc
```
