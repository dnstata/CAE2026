# Carregando o pacote necessário:
require(survival)


# Carregando dados:
tempos.75=c(8,109,219,423,514,973,1013,1100,1246,1503,1607,rep(1680,29))
tempos.95=c(8,67,155,160,161,196,204,308,518,523,551,560,639,734,802,855,943,1236,1266,1287,1307,
1351,1379,1521,1607,1639,1653,rep(1680,13))
tempos.115=c(1,11,11,12,24,32,59,63,75,81,84,86,130,133,135,142,142,195,202,202,215,224,224,239,
299,352,362,365,382,390,471,570,584,600,730,850,925,958,1413,1680)
tempos=c(tempos.75,tempos.95,tempos.115)
cens.75=c(rep(1,11),rep(0,29))
cens.95=c(rep(1,27),rep(0,13))
cens.115=c(rep(1,39),rep(0,1))
cens=c(cens.75,cens.95,cens.115)
tau.75=rep(75,40)
tau.95=rep(95,40)
tau.115=rep(115,40)
tau=c(tau.75,tau.95,tau.115)
x.75=1000/(273.15+tau.75)
x.95=1000/(273.15+tau.95)
x.115=1000/(273.15+tau.115)
x=c(x.75,x.95,x.115)


# Gráfico tempo de falha (h) versus níveis de temperatura (°C):
plot(tau, tempos, xlab="Temperatura (°C)", ylab="Tempo de falha (h)")


# Estimador K-M:
ekm=survfit(Surv(tempos,cens)~tau)
summary(ekm)


plot(ekm, main="Kaplan-Meier", xlab="Tempos (h)", ylab="R(t) estimada", mark.time=TRUE, col=c(1,2,4))
legend("bottomleft", c("75 °C","95 °C","115 °C"), col=c(1,2,4), bty="n", lty=1)


# Gráficos de linearização das funções de confiabilidade:
ekm.75=survfit(Surv(tempos.75,cens.75)~1)
t1=ekm.75$time
rt1=ekm.75$surv
invrt1=qnorm(rt1)


ekm.95=survfit(Surv(tempos.95,cens.95)~1)
t2=ekm.95$time
rt2=ekm.95$surv
invrt2=qnorm(rt2)


ekm.115=survfit(Surv(tempos.115,cens.115)~1)
t3=ekm.115$time
rt3=ekm.115$surv
invrt3=qnorm(rt3)


par(mfrow=c(1,4))
plot(t1, -log(rt1), pch=16, xlab="t", ylab="-log(R(t))", main="exponencial", xlim=c(0,1700), ylim=c(0,4), col=1)
lines(t2, -log(rt2), pch=16, type="p", col=2, xlim=c(0,1700), ylim=c(0,4))
lines(t3, -log(rt3), pch=16, type="p", col=4, xlim=c(0,1700), ylim=c(0,4))
legend("topleft", c("75 °C","95 °C","115 °C"), col=c(1,2,4), bty="n", lty=1)


plot(log(t1), log(-log(rt1)), pch=16, xlab="t", ylab="log(-log(R(t)))", main="Weibull",
xlim=c(0,8), ylim=c(-4,2), col=1)
lines(log(t2), log(-log(rt2)), pch=16, type="p", col=2, xlim=c(0,8), ylim=c(-4,2))
lines(log(t3), log(-log(rt3)), pch=16, type="p", col=4, xlim=c(0,8), ylim=c(-4,2))
legend("topleft", c("75 °C","95 °C","115 °C"), col=c(1,2,4), bty="n", lty=1)


plot(log(t1), invrt1, pch=16, xlab="log(t)", ylab=expression(Phiˆ -1*(R(t))), main="log-normal",
xlim=c(0,8), ylim=c(-2,2), col=1)
lines(log(t2), invrt2, pch=16, type="p", col=2, xlim=c(0,8), ylim=c(-2,2))
lines(log(t3), invrt3, pch=16, type="p", col=4, xlim=c(0,8), ylim=c(-2,2))
legend("bottomleft", c("75 °C","95 °C","115 °C"), col=c(1,2,4), bty="n", lty=1)


plot(log(t1),log(1/rt1-1), pch=16, xlab="log(t)", ylab="log(1/R(t)-1)", main="Dhillon",
xlim=c(0,8), ylim=c(-2,2), col=1)
lines(log(t2), log(1/rt2-1), pch=16, type="p", col=2, xlim=c(0,8), ylim=c(-2,2))
lines(log(t3), log(1/rt3-1), pch=16, type="p", col=4, xlim=c(0,8), ylim=c(-2,2))
legend("topleft", c("75 °C","95 °C","115 °C"), col=c(1,2,4), bty="n", lty=1)


# Ajuste dos modelos de regressão para todos os níveis de estresse:
mod.e=survreg(Surv(tempos,cens)~x, dist="exponential")
summary(mod.e)
AIC(mod.e)


mod.w=survreg(Surv(tempos,cens)~x, dist="weibull")
summary(mod.w)
AIC(mod.w)


gama.hat=1/mod.w$scale
gama.hat


mod.ln=survreg(Surv(tempos,cens)~x, dist="lognormal")
summary(mod.ln)
AIC(mod.ln)


mod.d=survreg(Surv(tempos, cens)~x, dist="loglogistic") 
summary(mod.d)
AIC(mod.d)


# Verificação da adequação dos modelos utilizando os resíduos de Cox-Snell:
par(mfrow=c(2,4))


b0.hat=mod.e$coefficients[2]
res=tempos*exp(-b0.hat-b1.hat*x)
ekm=survfit(Surv(res,cens)~1, type="kaplan-meier")
plot(ekm, conf.int=F, lty=c(1,1), xlab="e", ylab="R(e): estimada", main="exponencial")
res=sort(res); exp1=exp(-res); lines(res,exp1,lty=3)
legend("topright", lty=c(1,3), c("Kaplan-Meier","Exp(1)"), lwd=1, bty="n", cex=0.8)
rt=ekmtime; rexp1=exp(-t)
plot(rt, rexp1, xlab="R(e): Kaplan-Meier", ylab="R(e): Exp(1)", pch=16, main="exponencial")


b0.hat=mod.w$coefficients[2]; gama.hat=1/mod.w$scale
res=(tempos*exp(-b0.hat-b1.hat*x))**gama.hat
ekm=survfit(Surv(res,cens)~1, type="kaplan-meier")
plot(ekm, conf.int=F, lty=c(1,1), xlab="e", ylab="R(e): estimada", main="Weibull")
res=sort(res); exp1=exp(-res); lines(res,exp1,lty=3)
legend("topright", lty=c(1,3), c("Kaplan-Meier","Exp(1)"), lwd=1, bty="n", cex=0.8)
rt=ekm$time; rexp1=exp(-t)
plot(rt, rexp1, xlab="R(e): Kaplan-Meier", ylab="R(e): Exp(1)", pch=16, main="Weibull")


b0.hat=mod.ln$coefficients[2]; scale.hat=mod.ln$scale
res=-log(1-pnorm((log(tempos)-b0.hat-b1.hat*x)/scale.hat))
ekm=survfit(Surv(res,cens)~1, type="kaplan-meier")
plot(ekm, conf.int=F, lty=c(1,1), xlab="e", ylab="R(e): estimada", main="log-normal")
res=sort(res); exp1=exp(-res); lines(res,exp1,lty=3)
legend("topright", lty=c(1,3), c("Kaplan-Meier","Exp(1)"), lwd=1, bty="n", cex=0.8)
rt=ekmtime; rexp1=exp(-t)
plot(rt, rexp1, xlab="R(e): Kaplan-Meier", ylab="R(e): Exp(1)", pch=16, main="log-normal")


b0.hat=mod.d$coefficients[1]; b1.hat=mod.d$coefficients[2]; gama.hat=1/mod.d$scale
res=log(1 + exp(gama.hat*(log(tempos)-b0.hat-b1.hat*x)))
ekm=survfit(Surv(res,cens)~1, type="kaplan-meier")
plot(ekm, conf.int=F, lty=c(1,1), xlab="e", ylab="R(e): estimada", main="Dhillon")
res=sort(res); exp1=exp(-res); lines(res,exp1,lty=3)
legend("topright", lty=c(1,3), c("Kaplan-Meier","Exp(1)"), lwd=1, bty="n", cex=0.8)
rt=ekm$surv; t=ekm$time; rexp1=exp(-t)
plot(rt, rexp1, xlab="R(e): Kaplan-Meier", ylab="R(e): Exp(1)", pch=16, main="Dhillon")


# Extrapolação dos resultados para as condições de uso (55 °C):
# Considerando ainda o uso real de 5 horas/dia:


# MTTF:
b0.hat=mod.e$coefficients[2]
exp(b0.hat+b1.hat*1000/(273.16+55)) / 5


# Tempo mediano (t50):
-exp(b0.hat+b1.hat*1000/(273.16+55))*log(1-0.5) / 5


# Percentil 5% (t5):
-exp(b0.hat+b1.hat*1000/(273.16+55))*log(1-0.05) / 5


# 1-R(t) com t = 1 ano = 365*5 = 1825 horas (5 horas de uso, em média, por dia):
1-exp(-1825/exp(b0.hat+b1.hat*1000/(273.16+55)))


# Energia de ativação:
E=8.6171*10^(-5)*mod.e$coefficients[2]*1000
E