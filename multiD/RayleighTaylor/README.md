# Rayleigh-Taylor

このディレクトリには，Rayleigh-Taylor問題を解くサンプルコードと結果を可視化するためのファイルが入っている。

## 初期条件

RT不安定性をシミュレーションするためには，運動量保存式とエネルギー方程式に源項が必要である。重力加速度を $g$ とし $y$ 負方向に重力をかける( $g<0$ )と，

```math
\left(\frac{\partial \rho v_y}{\partial t}\right)_\mathrm{grav}
= \rho g
```

```math
\left(\frac{\partial E}{\partial t}\right)_\mathrm{grav}
= \rho g v_y
```

となる。


計算領域は $-1/4\le x\le 1/4,\quad -3/4\le y\le 3/4$ とする。
初期条件は，

```math
    \rho(x,y) = 
    \begin{cases}
        2 & \mathrm{for}~y\ge 0 \\
        1 & \mathrm{for}~y<0
    \end{cases}
```

```math
P(x,y) = P_0 + \rho(x,y) g y
```

```math
    v_y(x,y) = A \left\{1+\cos\left(2\pi \frac{y}{L_y}\right) \right\}
    \left\{-\cos\left(2\pi \frac{x}{L_x}\right)\right\} ,\quad v_x(x,y) = v_z(x,y) = 0
```

```math
    B_x(x,y) = B_0,\quad B_y(x,y) = 0,\quad B_z(x,y)=0
```

とする。重力加速度を $g=-0.1$ とし，境界面での圧力を $P_0=2.5$ とする。 $y$ 方向に一様な初期磁場を置く。摂動は初期不連続面を波数 $2\pi/L_x$ の波で揺らがせるように， $v_y$ に入れる。 $A=2.5\times 10^{-3}$ は初期摂動の大きさを表す。

$x$ 方向には周期境界条件を課す。$y$ 方向の境界条件は，揺らぎが0の場合に平衡状態 $\partial P/\partial y=\rho g$ を維持できるように設定する。サンプルコードでは以下の境界条件を設定している。

- $\rho$ と $v_x$ ， $v_z$ ， $\boldsymbol{B}$ , $\psi$ は勾配0境界にする。
- $v_y$ に対しては反射境界を課す。
- $P$ については，与えられた密度分布における平衡分布
- $\partial P/\partial y=\rho g$ を数値積分して代入する。


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

