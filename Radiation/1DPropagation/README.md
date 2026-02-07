# Radiation / 1DPropagation

このディレクトリには、輻射輸送の **1D 伝播（Propagation）** のサンプル計算を走らせて、出力を可視化するための一式が入っています。
Sn法、M1-closure法、Flux-limited diffusion法の3つの手法の比較ができます。

## コンパイル
解析サーバの場合、以下のような `Makefile`でコンパイルできます。
<pre>
##########################################
# Makefile for 1Dpropataion
##########################################

##########################################
# Programs
##########################################

exem1 = m1.x
exesn = sn.x
exefld = fld.x

######################	
# complile options
######################
# intel
#fc=ifort -diag-disable=10448
#fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

# GNU 
fc=gfortran
fopt=-g -fbacktrace -O2
#fopt= -O0 -g -fbacktrace -fcheck=all

all: ${exem1} ${exesn} ${exefld}

.PHONY: all clean allclean

#################
# simulation
#################

${exem1}: main_m1.f90
	mkdir -p m1mod
	${fc} ${fopt} $< -o $@
	rm *.mod

${exesn}: main_sn.f90
	${fc} ${fopt} $< -o $@
	rm *.mod

${exefld}: main_fld.f90
	${fc} ${fopt} $< -o $@
	rm *.mod

#################
# clean up
#################
clean:
	rm -f ${exem1} ${exesn} ${exefld} *.o *.mod *~ 
</pre>

## 実行
以下の実行と解析ではM1-closure法について説明します。Snの場合はm1をsnに、FLDの場合はm1をfldに読み替えてください。

	./m1.x
	
この実行で`m1/`というディレクトリが作成され、`m1/snap?????.dat`というファイルができます。バイナリで出力したい時は`main_m1.f90`の`flag_binary`を`.true.`にしてください。

## 可視化
簡単にリアルタイム解析をしたいときは`RealTimeAnim.plt`を実行します。
     
	 gnuplot RealTimeAnim.plt
	 
図をpng形式で保存したい場合は以下を行なってください。
     
	 gnuplot -e num=10 MakePngFile.plt
	 
このpngをAnimationにしたい場合以下を実行してください。
     
	 ./MakeAnim.plt m1 1 90
	 
