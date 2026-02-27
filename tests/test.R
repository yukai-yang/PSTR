# delete this file when the package is ok.
library(PSTR)

pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
               tvars=c('vala'), im=1, iT=14)

pstr

LinTest(pstr) 
print(pstr, mode="tests")

EstPSTR(use=pstr,im=1,iq=1,useDelta=T,par=c(-.57,-1.2), vLower=4, vUpper=4)
print(pstr, mode="estimates")

min(pstr$vg); max(pstr$vg)

EstPSTR(use=pstr,im=1,iq=1,useDelta=T,par=c(-0.462,0), method="CG")
print(pstr, mode="estimates")




ret <- plot_target(obj = pstr, iq = 1, basedon = c(1, 2),
                   from = c(-0.6, -4), to = c(-0.55, -2),
                   length.out = c(20, 20))

ret


pstr0 = pstr$clone()
EstPSTR(use=pstr0)
print(pstr0, mode="estimates")
