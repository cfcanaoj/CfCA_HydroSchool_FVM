# Dai-Woodward衝撃波管問題1

このディレクトリには，[Dai & Woodward ()]()が提案した衝撃波管問題を解くサンプルコードと結果を可視化するためのファイルが入っている。

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
1.08 \\
1.2 \\
0.01 \\
0.5 \\
0.95 \\
2/\sqrt{4\pi} \\
3.6/\sqrt{4\pi} \\
2/\sqrt{4\pi} \\
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
1.0 \\
0 \\
0 \\
0 \\
1.0 \\
2/\sqrt{4\pi} \\
4/\sqrt{4\pi} \\
2/\sqrt{4\pi} \\
\end{pmatrix},
$$

厳密解は以下の図のようになる。
<p align="center">
  <img src="./figures/daiwoodward1_exact.png" width="600">
</p>



## ディレクトリ内の構造

ファイル名など     | 説明  
------------------|----------
 main.f90     | サンプルコード
 ans/main.f90     | 完成品(実習のために，完成品のバックアップ)
 daiwoodward1_exact.dat | 衝撃波管問題の厳密解
 figures/daiwoodward1_exact.png | 衝撃波管問題の厳密解を示した図
 可視化スクリプト   | MakeCompare.py


## 可視化

### スナップショットの図示

**MakeCompare.py**

pythonを使って，複数の出力ディレクトリ内の特定snapXXXXX.datを比較するためのスクリプト

Usage: python MakeAnime.py [step_s] [step_e] [dir1] [dir2] ... 

lax/snap00040.datとhll/snap00040.datを比較する場合は，以下を実行する。
```bash
python3 MakeAnime.py 40 lax hdc
```
