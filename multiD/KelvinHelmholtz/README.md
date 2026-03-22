# Kelvin-Helmholtz不安定性

このディレクトリには，Kelvin-Helmholtz問題を解くサンプルコードと結果を可視化するためのファイルが入っている。

## 初期条件

計算領域を$-1/2\le x\le 1/2, -1\le y\le 1$とする。周期境界条件を$y$方向にも適用するために，接触面を$y=\pm 1/2$の2カ所に用意する。$|y|\le 1/2$の物理量を下付き添え字1を付けて表し，それ以外の領域の物理量を下付き添え字2を付けて表す。物理量を滑らかに変化させるため，以下のように$\tanh$をつかう。


$$
\boldsymbol{Q}(x,y) = \frac{\boldsymbol{Q}_1 - \boldsymbol{ Q}_2}{2}
\left[ 
\tanh\left( \frac{y+1/2}{h} \right) 
-\tanh\left( \frac{y-1/2}{h} \right) 
\right]
+\boldsymbol{Q}_2
$$

流体1と流体2の混合具合を調べるために，流体の速度で移流するスカラー場$C$を導入する。

$$
\frac{DC}{Dt}=
    \frac{\partial C}{\partial t}
     + \bm{v}\cdot\bm{\nabla} C = 0
     \;\;\;\Longrightarrow\;\;\;
    \frac{\partial \rho C}{\partial t}
     + \bm{\nabla} \cdot(\rho C\bm{v}) = 0
$$

$$
\boldsymbol{Q}_1 =
\begin{bmatrix}
\rho_1 \\
v_{x1} \\
v_{y1} \\
P_1 \\
C_1
\end{bmatrix}
=
\begin{bmatrix}
1 \\
1 \\
0 \\
1 \\
1
\end{bmatrix},
\qquad
\boldsymbol{Q}_2 =
\begin{bmatrix}
\rho_2 \\
v_{x2} \\
v_{y2} \\
P_2 \\
C_2
\end{bmatrix}
=
\begin{bmatrix}
1 \\
-1 \\
0 \\
1 \\
0
\end{bmatrix}
$$

流体1と流体2ともに音速は$\sqrt{\gamma}$である。


## ディレクトリ内の構造

ファイル名など     | 説明  
------------------|----------
 main_HDC.f90     | HDC法を実装したサンプルコード
 main_CT.f90      | CT法を実装したサンプルコード
 可視化スクリプト   | MakeAnime.py, MakeAnima_vars.py, MakeCompare.py, MakeCompare_vars.py


## 可視化

### 動画作成

**MakeAnime.py**

pythonを使って，特定の範囲のスナップショットから動画を作成するスクリプト。

```bash
> python3 MakeAnime.py -h
usage: python3 MakeAnime.py [ascii|binary] [step_s] [step_e] [dir1] [dir2] ... [-h] [--vmin VMIN] [--vmax VMAX] [--interval INTERVAL]

Create a comparison movie from multiple simulation outputs.

positional arguments:
  {ascii,binary}  input file format
  step_s          starting step number
  step_e          ending step number
  dirnames        one or more directories containing snapshot files

options:
  -h, --help      show this help message and exit

Example:
  python3 MakeAnime.py ascii 0 20 ct hdc
```

**MakeAnime_vars.py**

MakeAnime.pyに，他の物理量を表示できるように改良したスクリプト。

```bash
>python3 MakeAnime_vars.py -h
usage: python3 MakeAnime.py [ascii|binary] [step_s] [step_e] [varname] [linear|log] [dir1] [dir2] ... [-h] [--vmin VMIN] [--vmax VMAX] [--interval INTERVAL]

Create a comparison movie from multiple simulation outputs.

positional arguments:
  {ascii,binary}       input file format
  step_s               starting step number
  step_e               ending step number
  varname              variable to plot (rho, vx, vy, vz, P, Bx, By, Bz, Bpre, Ekin, beta)
  {linear,log}         color scale type
  dirnames             one or more directories containing snapshot files

options:
  -h, --help           show this help message and exit
  --vmin VMIN          (optional) manual minimum value for the color scale
  --vmax VMAX          (optional) manual maximum value for the color scale
  --interval INTERVAL  (optional) frame interval in milliseconds

Example:
  python3 MakeAnime.py ascii 0 20 rho linear ct hdc
  python3 MakeAnime.py binary 10 50 beta log ct hdc --vmin 1e-2 --vmax 1e2
```



### スナップショットの図示

**MakeCompare.py**

pythonを使って，複数の出力ディレクトリ内の特定snapXXXXX.datを比較するためのスクリプト

```bash
> python3 MakeCompare.py -h
usage: python3 MakeCompare.py [ascii|binary] [step] [dir1] [dir2] ... [-h]

Create a comparison shapshot from multiple simulation outputs.

positional arguments:
  {ascii,binary}  input file format
  step            step number
  dirnames        one or more directories containing snapshot files

options:
  -h, --help      show this help message and exit

Example:
  python3 MakeCompare.py ascii 20 ct hdc
```

**MakeCompare_vars.py**

MakeComapre.pyに，他の物理量を表示できるように改良したスクリプト。

```bash
> python3 MakeComapre.py -h
usage: python3 MakeCompare.py [ascii|binary] [step] [varname] [linear|log] [dir1] [dir2] ... [-h] [--vmin VMIN] [--vmax VMAX]

Create a comparison shapshot from multiple simulation outputs.

positional arguments:
  {ascii,binary}  input file format
  step            step number
  varname         variable to plot (rho, vx, vy, vz, P, Bx, By, Bz, Bpre, Ekin, beta)
  {linear,log}    color scale type
  dirnames        one or more directories containing snapshot files

options:
  -h, --help      show this help message and exit
  --vmin VMIN     (optional) manual minimum value for the color scale
  --vmax VMAX     (optional) manual maximum value for the color scale

Example:
  python3 MakeCompare.py ascii 20 rho linear ct hdc
  python3 MakeCompare.py binary 50 beta log ct hdc --vmin 1e-2 --vmax 1e2
```

