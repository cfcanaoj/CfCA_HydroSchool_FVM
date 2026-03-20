# 有限振幅円偏光Alfven波の伝播

このディレクトリには，有限振幅円偏光Alfven波の伝播を解くサンプルコードと結果を可視化するためのファイルが入っている。

## 初期条件

非摂動状態として，静止した一様ガスを考える($\rho=\rho_0$, $v=0$, $P=P_0$)。一様な磁場$\bm{B}=(B_0,0,0)$を考える。
以下のように$+x$方向に伝播するAlfven波を解く。

Alfv\'en波の非圧縮性は，無限小の振幅でのみ成り立ち，有限振幅では摂動の2次の項に相当する磁気圧$\delta \bm{B}_\perp^2$の変動により， 縦波が励起され得るが，円偏光したAlfv\'en波では磁気圧が空間的に厳密に一定となり，有限振幅でも厳密解になる。
右($+x$)方向に伝播する円偏光したAlfv\'en波解は，$A$を磁場の振幅とすると，円偏光したAlfv\'en波の磁場は，

$$
\bm{B} = 
\begin{pmatrix}
B_{x} \\
B_{y} \\
B_{z} 
\end{pmatrix}
=\begin{pmatrix}
B_0 \\
AB_0\sin(kx - \omega t) \\
AB_0\cos(kx - \omega t) \\
\end{pmatrix},
$$

と表せる。ここで，$\omega/k = c_\mathrm{A} = B_0/\sqrt{\rho_0}$

速度は，　
$$
\bm{v} = 
\begin{pmatrix}
v_{x} \\
v_{y} \\
v_{z} 
\end{pmatrix}
=\begin{pmatrix}
0 \\
-Ac_\mathrm{A}\sin(kx - \omega t) \\
-Ac_\mathrm{A}\cos(kx - \omega t) \\
\end{pmatrix},
$$
となる。

$\rho_0=1, P_0 = 10, B_0 = 1$とする。


## ディレクトリ内の構造

ファイル名など     | 説明  
------------------|----------
 ans/main.f90     | 完成品 (流体学校の実習のため，AlfvenWave/直下にはサンプルコードがない。)
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
