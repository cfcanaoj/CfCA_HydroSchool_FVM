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

dtout = 0.005

set yrange [-1.1e-5:1.1e-5]
set xrange [-0.5:0.5]
set grid

do for [i = ibeg:ifin ] {
  inpfile=sprintf("%s/snap%05d.dat",model,i)

  # extract time from the first line
  command=sprintf("awk 'NR==1{print($2)}' %s", inpfile)
  time=system(command)

  set title sprintf("time = %s",time) 
  plot inpfile u 1:($2-1) ti model w lp ,\
        1e-5*sin(2.0*pi*(x-time)) ti "exact solution" w l 
  pause 0.5
}
