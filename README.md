# CfCA Hydro School FVM (FY2025)

国立天文台 CfCA（Center for Computational Astrophysics）流体学校（FY2025）向けの、**有限体積法（Finite Volume Method; FVM）** による流体・MHD 計算の実習用コードとテスト問題集です。  
「動く完成品を配布する」だけでなく、講義内容（CFL条件、保存則、リーマンソルバー、MUSCL/TVD、境界条件、多次元効果など）を **手元で再現・比較できる** ことを目的にしています。

---

## 収録内容（ディレクトリ）

- `1D/` : 1次元テスト
  - `Sod`（衝撃波管）
  - `AcousticWave`（音波）
  - `BrioWu`（MHD衝撃波管）
  - `DaiWoodward1`, `DaiWoodward2`
  - `AlfvenWave`（アルヴェン波）
- `multiD/` : 多次元テスト
  - `BlastWave`
  - `AlfvenWave`
  - `ShockTubeRot`（回転座標/斜め衝撃波チューブ等）
  - `KelvinHelmholtz`（KH不安定）
  - `RayleighTaylor`（RT不安定）
  - `DecayingTurbulence`（減衰乱流）
- `Radiation/` : 輻射輸送
- `vis3D/` : 可視化補助（3D/多次元データ向け）
- `documents/` : 講義・実習用の資料

（上の一覧はリポジトリ直下の表示に合わせています。）

---

## 想定ユーザー

- 大学院生〜研究者で、流体/MHD の数値計算を「使う」だけでなく「中身を理解したい」人
- 近似リーマン解法、制限関数、時間積分、CFL、保存性、拡散/分散誤差の影響を、実際にパラメータを振って確認したい人

---

## 必要環境

最小構成の例：

- Fortran コンパイラ（例：`gfortran` / `ifort` / `nvfortran` など）
- `make`
- 解析・作図：Python 3（`numpy`, `matplotlib` 等）または gnuplot  
  ※どの可視化方法が同梱されているかは各ケースのディレクトリを参照してください。


