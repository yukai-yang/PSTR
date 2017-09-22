library(tidyverse)

dat = dat[,1]
length(dat)
dat = cbind(dat,lag(dat),lag(dat,2),lag(dat,3),lag(dat,4),lag(dat,5),lag(dat,6),lag(dat,7),lag(dat,8),lag(dat,9),lag(dat,10))
dat = dat[11:nrow(dat),]
dim(dat); head(dat)
colnames(dat) = paste0("spot_",0:10)

sunspot = as_tibble(dat)

ggplot(sunspot,aes(y=spot_0,x=1:nrow(sunspot))) + geom_line()

pstr = NewPSTR(sunspot, dep='spot_0', im=1, indep=2:11, tvars=2:7, iT=270)
LinTest(use=pstr)

pstr = EstPSTR(use=pstr, im=1, iq=2, par=c(5.5,7.9), method='CG')
print(pstr,"estimates")

ggplot(tibble(vg=pstr$vg,vq=pstr$mQ[,pstr$iq]),aes(y=vg,x=vq)) +
  labs(x=2) + scale_x_log10() +
  geom_point()


plot_transition(pstr,title=3, subtitle=4, caption=5,logx=T)
plot_transition(pstr,color="red",size=2)
