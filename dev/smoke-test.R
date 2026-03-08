devtools::load_all()
data("Hansen99", package = "PSTR")

pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
               tvars=c('vala'), im=1, iT=14)

pstr

LinTest(pstr)

pstr$WCB_LinTest(iB = 20, parallel = FALSE)

old_opt <- options(PSTR.future.globals.maxSize = 4 * 1024^3)
pstr$WCB_LinTest(iB = 20, parallel = TRUE, cpus = 2)
options(old_opt)

print(pstr, mode="tests")

EstPSTR(use=pstr,im=1,iq=1,useDelta=T,par=c(-.57,-1.2), vLower=4, vUpper=4)
print(pstr, mode="estimates")

min(pstr$vg); max(pstr$vg)

EstPSTR(use=pstr,im=1,iq=1,useDelta=T,par=c(-0.462,0), method="CG")
print(pstr, mode="estimates")

WCB_TVTest(use=pstr,iB=10,parallel=F)
WCB_HETest(use=pstr,vq=as.matrix(Hansen99[,'vala'])[,1],iB=10,parallel=F)

old_opt <- options(PSTR.future.globals.maxSize = 4 * 1024^3)
## wild bootstrap time-varying evaluation test 
WCB_TVTest(use=pstr,iB=20,parallel=TRUE,cpus=2)
## wild bootstrap heterogeneity evaluation test
WCB_HETest(use=pstr,vq=as.matrix(Hansen99[,'vala'])[,1],iB=20,parallel=TRUE,cpus=2)
options(old_opt)

print(pstr, mode="evaluation")


ret <- plot_target(obj = pstr, iq = 1, basedon = c(1, 2),
                   from = c(-0.6, -4), to = c(-0.55, -2),
                   length.out = c(20, 20))

ret


pstr0 = pstr$clone()
EstPSTR(use=pstr0)
print(pstr0, mode="estimates")

#######

# Hansen
pstr = NewPSTR(Hansen99, dep = "inva",
  indep = c("vala", "debta", "cfa", "sales"),
  indep_k = c("cfa"),
  tvars = c('vala','debta','cfa','sales'),
  im = 1, iT = 14 )


pstr = NewPSTR( Hansen99, dep = 'inva',
  indep = c('vala', 'cfa'),
  indep_k = c('vala', 'cfa', 'sales'),
  tvars = c('vala','debta','cfa','sales'),
  im = 1, iT = 14)
pstr


LinTest(pstr) 
print(pstr, mode="tests")

old_max <- getOption("PSTR.future.globals.maxSize")
options(PSTR.future.globals.maxSize = 4 * 1024^3)
pstr$WCB_LinTest(iB = 1000, parallel = TRUE, cpus = 6)
options(PSTR.future.globals.maxSize = old_max)

EstPSTR(use=pstr,im=1,iq=1,useDelta=T,par=c(1,0), method="CG")
print(pstr, mode="estimates")

log(pstr$gamma); pstr$c; pstr$s2

EstPSTR(use=pstr,im=1,iq=1,useDelta=F,par=c(pstr$gamma+0.2,pstr$c-0.2), method="CG")
print(pstr, mode="estimates")
ret <- plot_target(pstr, iq = 1, from = c(log(pstr$gamma-0.2), pstr$c-0.2),
                   to = c(log(pstr$gamma+0.2), pstr$c+0.2),
                   length.out = c(30, 30))
ret

## finally
EstPSTR(use=pstr,im=1,iq=1,useDelta=T,par=c(2, 0), method="CG")
print(pstr, mode="estimates")

plot_transition(pstr)

#
pcoef <- plot_coefficients(pstr, vars = 1:4)
pcoef[[1]]

EvalTest(use=pstr, vq = as.matrix(Hansen99[,'vala'])[,1])
print(pstr, mode = "evaluation")

iB = 400
old_max <- getOption("PSTR.future.globals.maxSize")
options(PSTR.future.globals.maxSize = 4 * 1024^3)  # 4 GB, if needed
pstr <- WCB_TVTest(pstr, iB = iB, parallel = TRUE, cpus = 6)
pstr <- WCB_HETest(pstr, vq = as.matrix(Hansen99[, "vala"])[, 1],
  iB = iB, parallel = TRUE, cpus = 6)
# Reset to the previous option value.
options(PSTR.future.globals.maxSize = old_max)

print(pstr, mode = "evaluation")