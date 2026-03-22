# Blast waveテスト

このディレクトリには，Blast wave問題を解くサンプルコードと結果を可視化するためのファイルが入っている。

## 初期条件

Blast waveの問題設定を以下に示す。
一様磁場を持つ，一様で静止したガスの中心の狭い領域に，高圧のガスを置く。
計算を開始すると，高圧ガスが外に向かって膨張し，衝撃波が外に伝播する。

```math
    \rho=1, \quad\boldsymbol{v}=\boldsymbol{0},\quad\boldsymbol{B}=(B_0/\sqrt{2},B_0/\sqrt{2},0)
```

```math
P =
\begin{cases}
10   & \text{if } |\boldsymbol{r}| < 0.1 \\
0.1  & \text{if } |\boldsymbol{r}| \ge 0.1
\end{cases}
```


ここで，初期磁場強度は $B_0=10$ にする。

## ディレクトリ内の構造

ファイル名など     | 説明  
------------------|----------
 main_HDC.f90     | HDC法を実装したサンプルコード
 main_CT.f90      | CT法を実装したサンプルコード
 可視化スクリプト   | MakeAnime.py, MakeAnime_vars.py, MakeCompare.py, MakeCompare_vars.py


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

