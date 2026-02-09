# Radiation / 2D CrossBeams

このディレクトリには、輻射輸送の **二次元ビーム交差 （2D cross beams test）** のサンプル計算を走らせ、出力を可視化するための一式が入っています。  
**S\_N 法（SN） / M1-closure（M1） ** の 2手法を同一設定で比較できます。


## コンパイル
解析サーバの場合、以下でコンパイルできます。XD2000の場合は`Makefile`を編集してください。
```bash
make
```
生成される実行ファイルは`m1.x`（M1）,`sn.x`（SN）です。
以下のコマンド実行ファイル等を消せます。
```bash
make clean
```
※コンパイラ切り替え（ifort を使う等）は Makefile 内の fc を変更してください。

## 実行
以下では M1（`m1.x`） を例に説明します。
SN の場合は`m1`を`sn`に読み替えてください。
```bash
	./m1.x
```
この実行で`m1/`というディレクトリが作成され、`m1/snap?????.dat`というスナップショットが出力されます。
バイナリで出力したい時は`main_m1.f90`の`flag_binary`を`.true.`にしてください。

## 可視化 (gnuplot)
簡単にリアルタイム解析をしたいときは`RealTimeAnim.plt`を実行します。
```bash
gnuplot RealTimeAnim.plt
```
図をpng形式で保存したい場合は以下を行なってください。`num` は保存するスナップショット番号（例: `snap00010.dat` 相当）です。
```bash
gnuplot -e num=10 MakePngFile.plt
```	 
このpngをAnimationにしたい場合以下を実行してください。第1引数: 対象ディレクトリ　第2,3引数: 使うスナップショット番号の開始と終了
```bash
./MakeAnim.plt m1 1 90
```
