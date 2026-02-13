set terminal x11

if (ARGC < 3) {
   print "Usage: gnuplot -c RealtimeAnim.plt model ibeg ifin"
   print "Now set model = lax, ibeg = 1, ifin = 40"
   model="lax"
   ibeg=1
   ifin=40
} else {
   model = (ARG1)
   ibeg = int(ARG2)
   ifin = int(ARG3)
}
set xlabel "x"

set yrange [0.0:1.2]
set xrange [-0.5:0.5]

do for [i = ibeg:ifin ] {
  inpfile=sprintf("%s/snap%05d.dat",model,i)

  # extract time from the first line
  command=sprintf("awk 'NR==1{print($2)}' %s", inpfile)
  time=system(command)

  set title sprintf("time = %s",time) 
  plot inpfile u 1:2 ti "numerical" w lp ,\
      "sod_ana.dat" u ($1*time/0.2):2 ti "exact solution" w l
  pause 0.1
}
