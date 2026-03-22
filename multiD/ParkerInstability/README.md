# Parker不安定性

このディレクトリには，Rayleigh-Taylor問題を解くサンプルコードと結果を可視化するためのファイルが入っている。

## 初期条件

重力加速度は，中心面 $y=0$ で0となり，$|y|$ が増大し，スケールハイト $H_g$ を超えると，大きさが $g_0$ に漸近するような形として，以下の表式を採用する。

```math
    \boldsymbol{g} = \left(0, - g_0 \tanh(y/H_g)\right)
```

重力加速度の方向は中心面 $y=0$ を向いている。
この設定により， $y$ 方向境界を中心面付近で発展するParker不安定性から
十分離した位置に配置でき，境界の影響を最小限に抑えられる。

磁場は $x$ 方向を向き，磁気圧と圧力の比(プラズマベータ)が一定であるとする。

```math
    \boldsymbol{B} = (B_x(y),0,0),\quad \beta_0 = \frac{2P(y)}{B_x(y)^2} = \mathrm{const.}
```

鉛直方向の静水圧平衡の式は，

```math
\frac{d}{dy}\left( P + \frac{B_x^2}{2}\right) = \rho g_y\quad
\Longrightarrow\quad
  (1+\beta_0^{-1})  \frac{dP}{dy} = \rho g_y
```

となる。温度 $T\equiv P/\rho$ は， $|y|<y_0$ で低温ガス(温度 $T_\mathrm{L}$ )， $|y|>y_0$ で高温ガス( $T_\mathrm{H}$ )となるように設定し，両者をスケール $H_t$ で滑らかに接続する。

```math
    T = T_\mathrm{L} + \frac{1}{2} (T_\mathrm{H} - T_\mathrm{L})\left[
    \tanh\left(\frac{|y|-y_0}{H_t}\right)+1
    \right]
```

この式を使うと，静水圧平衡の式は，圧力に関する常微分方程式となるので，数値的に積分できる。以下の形式解が得られる。

```math
    P = P(y=0)
    \exp\left[
    \int_0^y \frac{g_y}{(1+\beta_0^{-1})T}dy
    \right]
```

$y=0$ における境界条件は， $\rho(y=0)=1$ である。 $y=0$ において音速が1となるように，低温ガスの温度は $T_\mathrm{L} = 1/\gamma$ とする。上空の温度 $T_\mathrm{H}$ は $25T_\mathrm{L}$ とし， $y_0=10$ において，スケール $H_t=0.5$ で低温から高温へ遷移する。プラズマベータは $\beta_0=1$ とする。重力加速度の漸近値は  $g_0=1.47$ ，重力のスケール長は $H_g=5$ とする。比熱比は $\gamma=1.05$ とする。

計算領域は $|x|\le 7.5\pi$ ， $|y|\le 15\pi$ とする。

上記の静水圧平衡状態に，水平方向の速度擾乱を加える。

```math
 \delta v_x = A\sin\left(\frac{2\pi x}{\lambda}\right)
    \frac{1}{2}\left\{\left[ \tanh\left(\frac{y+4}{0.5}\right)-\tanh\left(\frac{y+1}{0.5}\right)\right]
+\left[\tanh\left(\frac{y-4}{0.5}\right)-\tanh\left(\frac{y-1}{0.5}\right)\right]\right\}
```

速度擾乱は， $1<|y|<4$ 付近に集中するように設定されており，ある $x$ において， $y>0$ と $y<0$ で $\delta v_x$ の符号が反転し，中心面に対して奇関数となる。波長は $\lambda=7.5\pi$ とする。振幅 $A$ は $10^{-2}$ とする。

 $x$ 方向には周期境界条件を課し， $y$ 方向の境界条件は，初期値に固定するように置く。


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

