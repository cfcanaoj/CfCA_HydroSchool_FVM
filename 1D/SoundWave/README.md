# 音波の伝播

このディレクトリには，音波の伝播を解くサンプルコードと結果を可視化するためのファイルが入っている。

## 初期条件

非摂動状態として，静止した一様ガスを考える($\rho=\rho_0=1$, $v=0$, $P=P_0=1/\gamma$)。
以下のように$+x$方向に伝播する密度摂動を考える。
$$
\rho = \rho_0 + \delta \rho\sin(k(x - c_\mathrm{s}t)),
$$
ここで，$k=2\pi/L$は波数, $L$は解散領域のサイズ，$c_\mathrm{s}=\sqrt{\gamma P_0 /\rho_0}$は音速である。

上記の密度摂動と整合的な速度摂動と圧力摂動を求めよう。
線形化した連続の式
$$
\rho_0 \frac{\partial v}{\partial x} + \frac{\partial \rho}{\partial t} = 0
$$
から，以下のように速度摂動が導かれる。
$$\rho_0 \frac{\partial \delta v}{\partial x}
= c_\mathrm{s}k \delta \rho\cos (k(x-c_\mathrm{s}t)),\;\;\;
\Rightarrow\;\;\;
\delta v = c_\mathrm{s} \frac{\delta \rho}{\rho_0} \sin (k(x-c_\mathrm{s}t))
$$
断熱の式を線形化すると($d\ln (P\rho^{-\gamma})/dt = 0$), 以下のように圧力摂動を得る。
$$
\frac{\delta P}{P_0} = \gamma \frac{\delta \rho}{\rho_0}\sin(k(x - c_\mathrm{s}t)),
\;\;\;\Rightarrow\;\;\;
\delta P = c_\mathrm{s}^2 \rho_0 \sin(k(x-c_\mathrm{s}t))
$$

計算領域は$-0.5<x<0.5$とし，計算領域内に1波長入る波数を考える($k=2\pi$)。



## ディレクトリ内の構造

ファイル名         | 説明  
------------------|----------
 ans/main.f90     | 完成品
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

* MakeCompare.py

pythonを使って，複数の出力ディレクトリ内の特定snapXXXXX.datを比較するためのスクリプト

Usage: python MakeAnime.py [step_s] [step_e] [dir1] [dir2] ... 

lax/snap00040.datとhll/snap00040.datを比較する場合は，以下を実行する。
```bash
python3 MakeAnime.py 40 lax hdc
```
