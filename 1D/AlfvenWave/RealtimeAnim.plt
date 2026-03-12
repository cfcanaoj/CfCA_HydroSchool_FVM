set terminal x11

if (ARGC < 3) {
   print "Usage: gnuplot -c RealtimeAnim.plt ibeg ifin model"
   print "Now set ibeg = 1, ifin = 40, model = lax"
   ibeg=1
   ifin=40
   model="lax"
} else {
   ibeg = int(ARG1)
   ifin = int(ARG2)
   model = (ARG3)
}
set xlabel "x"

set yrange [-0.12:0.12]
set xrange [-0.5:0.5]

do for [i = ibeg:ifin ] {
  inpfile=sprintf("%s/snap%05d.dat",model,i)

  # extract time from the first line
  command=sprintf("awk 'NR==1{print($2)}' %s", inpfile)
  time=system(command)

  set title sprintf("time = %s",time) 
  plot inpfile u 1:8 ti model w lp ,\
       1e-1*sin(2.0*pi*(x-time)) ti "exact solution" w l 
  pause 0.1
}
