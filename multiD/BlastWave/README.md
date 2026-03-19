# Blast Wave Test

## Description

## Setup
The initial condition is 

$$ {\left(\matrix{\rho_{\rm L} \cr v_{\rm L} \cr P_{\rm L} } \right)}
= {\left(\matrix{1 \cr 0 \cr 1 } \right)},
{\left(\matrix{\rho_{\rm R} \cr v_{\rm R} \cr P_{\rm R} } \right)}
= {\left(\matrix{0.125 \cr 0 \cr 0.1 } \right)}
$$

## Extension

## Visualization

* RealtimeAnim.plt

`gnuplot -c RealtimeAnim.plt ibeg ifin dirname`

* MakeAnime.sh

`./MakeAnime.sh ibeg ifin dirname`


* MakeAnim.py

Usage: python MakeAnime.py [step_s] [step_e] [dir1] [dir2] ... 

* MakeCompare.py

usage: python3 MakeCompare.py [step] [dir1] [dir2] ...
