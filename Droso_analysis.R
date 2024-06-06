copudata <- read.csv("C:/Users/sumit/OneDrive/Documents/Drosophila manuscript/Droso_analysis_orig/Copulation duration data.csv")
copudata$grp <- paste(copudata$Female.age, copudata$Male.age, sep='.')

#group statistics
grp.mean  <- tapply(copudata$Copulation.duration.seconds, list(copudata$Female.age, copudata$Male.age), mean)
grp.sd  <- tapply(copudata$Copulation.duration.seconds, list(copudata$Female.age, copudata$Male.age), sd)

gmeans <- data.frame(means=as.vector(grp.mean), Female.age=rep(c(4,11),5), Male.age=rep(c(4,8,11,15,19), rep(2,5)))
gsd <- data.frame(sd=as.vector(grp.sd), Female.age=rep(c(4,11),5), Male.age=rep(c(4,8,11,15,19), rep(2,5)))
copudata2 <- merge(copudata, gmeans)
copudata3 <- merge(copudata2, gsd)

#Figure 1
library(lattice)
dotplot(factor(Male.age)~Copulation.duration.seconds|factor(Female.age, levels=sort(unique(copudata3$Female.age)), labels=paste('Female age: ', sort(unique(copudata3$Female.age)), ' days', sep='')), data=copudata3, layout=c(1,2), ylab='Age of male (days)', xlab='Duration of copulation (seconds)', panel=function(x, y, subscripts){
panel.dotplot(x, y, pch=1, cex=.6)
panel.arrows(copudata3$means[subscripts]-copudata3$sd[subscripts], as.numeric(y)+.2, copudata3$means[subscripts] + copudata3$sd[subscripts], as.numeric(y)+.2, code=3, angle=90, length=.05, col=3)
panel.points(copudata3$means[subscripts], as.numeric(y)+.2, pch=4, col=3)
})

#Figure 2

plot(as.vector(grp.sd)~as.vector(grp.mean), xlab='Mean', ylab='Standard deviation')
abline(lm(as.vector(grp.sd)~as.vector(grp.mean)), lty=2)
lines(lowess(as.vector(grp.sd)~as.vector(grp.mean), f=.8), col=2)
legend('topleft', c('linear regression', 'lowess'), lty=c(2,1), col=c(1,2), cex=.9, bty='n')

#Table 1

out1 <- lm(Copulation.duration.seconds~factor(Male.age) * factor(Female.age), data=copudata)
out3 <- glm(Copulation.duration.seconds~factor(Male.age) * factor(Female.age), data=copudata, family=Gamma)
out2 <- glm(Copulation.duration.seconds~factor(Male.age) * factor(Female.age), data=copudata, family=Gamma(link=log))
AIC(out1,out3,out2)
out4 <- glm(Copulation.duration.seconds~Male.age * factor(Female.age), data=copudata, family=Gamma(link=log))
AIC(out1,out3,out2)

anova(out2, test='Chisq')
summary(out2)
summary(out4)

Anova(out2,type=3)


# Figure 3

out2a <- glm(Copulation.duration.seconds~factor(Male.age) * factor(Female.age), data=copudata3, family=Gamma(link=log))
copudata3$gam.mean <- fitted(out2a)
myrate <- summary(out2a)$dispersion
copudata3$gam.shape <- 1/myrate
copudata3$scale <- copudata3$gam.mean*myrate

my.panel <- function(x, y, subscripts, col, lty, pch, group.number, ...) {
loc <- c(-.0005, -.0012)
panel.mathdensity(dmath=dgamma, args=list(shape=copudata3$gam.shape[1], scale=y), col=col, lty=lty, n=100)
#save old graphics settings and restore them on exit
opar <- trellis.par.get()
on.exit(trellis.par.set(opar))
trellis.par.set(list(box.rectangle=list(col=c('grey50', 'orangered2')[group.number], lty=lty), box.umbrella=list(col=c('grey50', 'orangered2')[group.number], lty=lty)))
panel.bwplot(x, rep(loc[group.number], length(x)), horizontal=T, col=col, box.width=.0005, pch='|')
panel.xyplot(x, jitter(rep(loc[group.number], length(x)), a=.00025), col=col, cex=.45, pch=pch)
}

xyplot(scale ~ Copulation.duration.seconds|factor(copudata3$Male.age, levels=sort(unique(copudata3$Male.age)), labels=paste('Male age =', sort(unique(copudata3$Male.age)), 'days')), data=copudata3, ylim=c(-.0017, 0.0085), ylab='Probability Density', xlab='Duration of copulation (seconds)', groups=factor(Female.age, levels=sort(unique(copudata3$Female.age))), col=c(1,'red3'), lty=c(1,2), pch=c(15,16), panel.groups='my.panel', panel="panel.superpose", key=list(x=.78, y=.84, corner=c(0,0), points=list(pch=c(15,16), col=c(1,'red3'), cex=.75), text=list(paste(sort(unique(copudata3$Female.age)), 'days'), cex=.85), lines=list(col=1:2, lty=1:2, size=3),  title='Female age', cex.title=.9), par.settings=list(plot.symbol = list(pch=1, col="white", cex=.01)))

# Figure 4

out2b <- glm(Copulation.duration.seconds ~ factor(Male.age):factor(Female.age)-1, data=copudata3, family=Gamma(link=log))

out.coef <- summary(out2b)$coefficients
out.vcov <- vcov(out2b)
out.mod <- data.frame(expand.grid(1:5, 1:2), est=out.coef[,1], se=out.coef[,2])
out.mod$low95 <- out.mod$est - qnorm(.975) * out.mod$se
out.mod$up95 <- out.mod$est + qnorm(.975) * out.mod$se

nor.func1 <- function(alpha, model, sig) 1 - pt(-qt(1-alpha/2, model$df.residual) * sum(sqrt(diag(sig))) / sqrt(c(1,-1) %*% sig %*%c (1,-1)), model$df.residual) - pt(qt(1-alpha/2, model$df.residual) * sum(sqrt(diag(sig))) / sqrt(c(1,-1) %*% sig %*% c(1,-1)), model$df.residual, lower.tail=F)
nor.func2 <- function(a,model,sigma) nor.func1(a, model, sigma)-.95
ci.func <- function(rowvals, glm.model, glm.vmat) {
n <- length(rowvals)
xvec1b <- numeric(n*(n-1)/2)
vmat <- glm.vmat[rowvals, rowvals]
ind <- 1
for(i in 1:(n-1)) {
for(j in (i+1):n){
sig <- vmat[c(i,j), c(i,j)]
#solve for alpha
xvec1b[ind] <- uniroot(function(x) nor.func2(x, glm.model, sig), c(.001,.999))$root
ind <- ind+1
}}
1-xvec1b
}

vals1 <- ci.func(1:5, out2b, out.vcov)
vals2 <- ci.func(6:10, out2b, out.vcov)
vals1
vals2

ci.val <- c(rep(mean(vals1),5), rep(mean(vals2),5))
out.mod$low50 <- out.mod$est + out.mod$se * qnorm((1-ci.val)/2)
out.mod$up50 <- out.mod$est + out.mod$se * qnorm(1-(1-ci.val)/2)

myprepanel.ci <- function(x,y,lx,ux,subscripts,...) {
list(xlim=range(x, ux[subscripts], lx[subscripts], finite=TRUE))
}

dotplot(factor(Var1, levels=1:5, labels=c(4,8,11,15,19)) ~ exp(est)|factor(Var2, levels=1:2, labels=paste('Female age = ', c(4,11), ' days', sep='')), data=out.mod, xlab='Mean copulation duration time (seconds)', ylab='Male age (days)', prepanel=myprepanel.ci, lx=exp(out.mod$low95), ux=exp(out.mod$up95), layout=c(2,1), panel=function(x, y, subscripts){
panel.dotplot(x, y, col=4, pch='+', cex=.6)
panel.segments(exp(out.mod$low95)[subscripts], y, exp(out.mod$up95)[subscripts], y, lwd=3, col='dodgerblue4', lineend=1)
panel.segments(exp(out.mod$low50)[subscripts], y, exp(out.mod$up50)[subscripts], y, lwd=6, col='dodgerblue1', lineend=1)
panel.points(x, y, col='white', pch=16, cex=1.1)
panel.points(x, y, col='dodgerblue4', pch=1, cex=1.2)
panel.abline(v=0, col=2, lty=2)
}, scales=list(x='free'), strip=strip.custom(par.strip.text = list(cex=0.9)))

# Table 2
out3 <- glm(formula = Copulation.duration.seconds ~ Male.age * factor(Female.age) + I(Male.age^2) * factor(Female.age), family = Gamma(link = log), data = copudata)
AIC(out2a, out3)
sapply(list(out2a, out3), logLik)

# Figure 5
out3a <- glm(formula = Copulation.duration.seconds ~ Male.age : factor(Female.age) + I(Male.age^2) : factor(Female.age) + factor(Female.age)-1, family = Gamma(link = log), data = copudata)
summary(out3a)

func.fem11 <- function(x) exp(coef(out3a)[2] + coef(out3a)[4]*x + coef(out3a)[6]*x^2)
func.fem4 <- function(x) exp(coef(out3a)[1] + coef(out3a)[3]*x + coef(out3a)[5]*x^2)
curve(func.fem4, from=3, to=20, xlab='Male age (days)', ylab='Mean copulation duration (secs)')
curve(func.fem11, add=T, col=2)
points(as.numeric(colnames(grp.mean)), grp.mean[1,], col=1)
points(as.numeric(colnames(grp.mean)), grp.mean[2,], col=2, pch=16)
legend('topleft', c('4 days','11 days'), col=c(1,2), lty=1, pch=c(1,16), cex=.9, title='Female age', bty='n')

#Table 3

eggdata <- read.csv("/Users/Sumit/Desktop/Droso Analyt from Jack/Egg data.csv")
eggdata$grp <- paste(eggdata$Female.age, eggdata$Male.age, sep='.')
eggdata$grp.f <- factor(eggdata$grp, level=c('4.4', '4.11', '4.19', '11.4', '11.11', '11.19'))



egg.mean <- tapply(eggdata$Eggs, eggdata$grp, mean)
egg.sd <- tapply(eggdata$Eggs, eggdata$grp, sd)
egg.n <- tapply(eggdata$Eggs, eggdata$grp, length)
egg.var <- tapply(eggdata$Eggs, eggdata$grp, var)

library(MASS)
all.nb <- lapply(unique(eggdata$grp), function(x) {
  cur.dat <- eggdata[eggdata$grp==x,]
  fitdistr(cur.dat$Eggs, "negative binomial")})

all.pois <- lapply(unique(eggdata$grp), function(x) {
  cur.dat <- eggdata[eggdata$grp==x,]
  fitdistr(cur.dat$Eggs, "Poisson")})

egg.stats <- data.frame(Female.age=rep(c(4,1), c(3,3)), Male.age=rep(c(4,11,19),2), mean=egg.mean, var=egg.var, ratio=egg.var/egg.mean, Pois.LL=sapply(all.pois,function(x) x$loglik), NB.LL=sapply(all.nb,function(x) x$loglik))
egg.stats

# Figure 6

plot(egg.var~egg.mean, xlab='Mean', ylab='Variance')
abline(lm(egg.var~egg.mean), lty=2)
mod2 <- lm(egg.var~egg.mean + I(egg.mean^2))
curve(coef(mod2)[1] +coef(mod2)[2]*x + coef(mod2)[3]*x^2, add=T, col=2, lty=2)
anova(mod2)

# Table 4

nb.stats <- data.frame(Female.age=rep(c(4,1), c(3,3)), Male.age=rep(c(4,11,19),2), t(sapply(1:6, function(x) all.nb[[x]]$estimate)))
nb.stats

# Figure 7

out.nb2 <- data.frame(t(sapply(1:6, function(x) all.nb[[x]]$estimate)), grp=unique(eggdata$grp))
eggdata2 <- merge(eggdata, out.nb2)

myplot <- xyplot(Eggs~mu|factor(Male.age, levels=c(4,11,19), labels=paste('Male age = ', c(4,11,19), sep='')) + factor(Female.age, levels=c(4,11), labels=paste('Female age = ', c(4,11), sep='')), data=eggdata2, xlim=c(-15,250), ylim=c(-.005,.02), layout=c(3,2), xlab='Number of eggs', ylab='Probability', panel=function(x, y, subscripts){
cur.mu <- eggdata2$mu[subscripts][1]
cur.size <- eggdata2$size[subscripts][1]
cur.y <- dnbinom(0:250, mu=cur.mu, size=cur.size)
panel.points(0:250, cur.y, type='h', col='grey80')
panel.bwplot(y, rep(-.0025, length(x)), horizontal=T,  box.width=.003, col=1, pch='|')
panel.xyplot(y, jitter(rep(-.0025, length(x)), a=.001), cex=.45, pch=16, col='grey40')
}, par.settings=list(box.rectangle = list(col="grey50", lty=1), box.umbrella = list(col="grey50", lty=1), plot.symbol = list(pch=1, col="white", cex=.01)))

library(latticeExtra)
useOuterStrips(myplot)

#Figure 8

boxplot(Eggs~grp.f, data=eggdata, outline=F, axes=F, ylab='Number of eggs', xlab='Male age')
points(jitter(as.numeric(eggdata$grp.f), a=.2), eggdata$Eggs, pch=16, cex=.7, col='seagreen')
axis(1, at=1:6, labels=rep(c(4,11,19),2))
axis(2)
box()
mtext(at=2, side=3, 'Female age = 4', line=.4, cex=.9)
mtext(at=5, side=3, 'Female age = 11', line=.4, cex=.9)
abline(v=3.5, col=2, lty=5)
points(1:6, tapply(eggdata$Eggs, eggdata$grp.f, mean), col=2, pch=8, cex=1.25)

# Table 5

library(gamlss)
out.g0 <- gamlss(Eggs~1, data=eggdata, family=NBI)
out.g0a <- gamlss(Eggs~1,sigma.formula=~grp.f, data=eggdata, family=NBI)
out.g1 <- gamlss(Eggs~factor(Male.age), data=eggdata, family=NBI)
out.g1a <- gamlss(Eggs~factor(Male.age), sigma.formula=~grp.f, data=eggdata, family=NBI)
out.g2 <- gamlss(Eggs~factor(Female.age), data=eggdata, family=NBI)
out.g2a <- gamlss(Eggs~factor(Female.age), sigma.formula=~grp.f, data=eggdata, family=NBI)
out.g3 <- gamlss(Eggs~factor(Female.age) + factor(Male.age), data=eggdata, family=NBI)
out.g3a <- gamlss(Eggs~factor(Female.age) + factor(Male.age), sigma.formula=~grp.f, data=eggdata, family=NBI)
out.g4 <- gamlss(Eggs~factor(Female.age) * factor(Male.age), data=eggdata, family=NBI)
out.g4a <- gamlss(Eggs~factor(Female.age) * factor(Male.age), sigma.formula=~grp.f, data=eggdata, family=NBI)

data.frame(k=sapply(list(out.g0, out.g1, out.g2, out.g3, out.g4, out.g0a, out.g1a, out.g2a, out.g3a, out.g4a), function(x) attr(logLik(x), 'df')), LL=round(sapply(list(out.g0, out.g1, out.g2, out.g3, out.g4, out.g0a, out.g1a, out.g2a, out.g3a, out.g4a), logLik), 1), AIC=round(sapply(list(out.g0, out.g1, out.g2, out.g3, out.g4, out.g0a, out.g1a, out.g2a, out.g3a, out.g4a), AIC), 1))

#test female age
LR.test(out.g0a, out.g2a)

#test female age controlling for male age
LR.test(out.g1a, out.g3a)

summary(out.g3a)

# Table 6

nb9a <- gamlss(Eggs~Copulation.duration.seconds * factor(Male.age) * factor(Female.age), data=eggdata, family=NBII)
nb8a <- gamlss(Eggs~(Copulation.duration.seconds + factor(Male.age) + factor(Female.age))^2, data=eggdata, family=NBII)
LR.test(nb8a, nb9a)
nb7a <- gamlss(Eggs~Copulation.duration.seconds*factor(Male.age) + Copulation.duration.seconds*factor(Female.age), data=eggdata, family=NBII)
LR.test(nb7a, nb8a)
nb6a <- gamlss(Eggs~Copulation.duration.seconds*factor(Male.age) + factor(Female.age), data=eggdata,family=NBII)
LR.test(nb6a, nb7a)
nb5a <- gamlss(Eggs~Copulation.duration.seconds + factor(Male.age) + factor(Female.age), data=eggdata, family=NBII)
LR.test(nb5a, nb6a)


sapply(list(nb9a,nb8a,nb7a,nb6a,nb5a), AIC)
sapply(list(nb9a,nb8a,nb7a,nb6a,nb5a), logLik)
sapply(list(nb9a,nb8a,nb7a,nb6a,nb5a), function(x) attr(logLik(x), 'df'))

summary(nb6a)

# Figure 9
min.cop <- tapply(eggdata$Copulation.duration.seconds, list(eggdata$Male.age, eggdata$Female.age), min)
max.cop <- tapply(eggdata$Copulation.duration.seconds, list(eggdata$Male.age, eggdata$Female.age), max)
eggdata$Male.age.f <- factor(eggdata$Male.age, levels=unique(eggdata$Male.age))

mypch <- c(15, 16, 17, 22, 21, 24)
set.seed(10)
plot(Eggs~jitter(Copulation.duration.seconds,a=5), data=eggdata, pch=mypch[eggdata$grp.f], cex=.8, col=c(1,2,4)[eggdata$Male.age.f], xlab='Copulation duration time (seconds)', ylab='Number of eggs', ylim=c(0,230))
nb.model <- function(x,Male,Female) exp(coef(nb6a)[1] + coef(nb6a)[2]*x + coef(nb6a)[3]*(Male==11) + coef(nb6a)[4]*(Male==19) + coef(nb6a)[5]*(Female==11) + coef(nb6a)[6]*x*(Male==11) + coef(nb6a)[7]*x*(Male==19))
curve(nb.model(x,4,4), add=T, from=min.cop[1,1], to=max.cop[1,1])
curve(nb.model(x,4,11), add=T, from=min.cop[1,2], to=max.cop[1,2], lty=2)
curve(nb.model(x,11,4), add=T, from=min.cop[2,1], to=max.cop[2,1], col=2)
curve(nb.model(x,11,11), add=T, from=min.cop[2,2], to=max.cop[2,2], col=2, lty=2)
curve(nb.model(x,19,4), add=T, from=min.cop[3,1], to=max.cop[3,1], col=4)
curve(nb.model(x,19,11), add=T, from=min.cop[3,2], to=max.cop[3,2], col=4,lty=2)
legend(80,235, c('4 day male','11 day male', '19 day male'), pch=mypch[1:3], col=c(1,2,4), lty=1, cex=.8, bty='n', pt.cex=1, seg.len=2.5, title=expression(bold('4 day female')))
legend(250,235, c('4 day male', '11 day male', '19 day male'), pch=mypch[4:6], col=c(1,2,4), lty=2, cex=.8, bty='n', pt.cex=1, seg.len=2.5, title=expression(bold('11 day female')))
 
#verify slope is zero for two male ages
nb6a.ex <- gamlss(Eggs~Copulation.duration.seconds:factor(Male.age) + factor(Male.age) + factor(Female.age)-1, data=eggdata, family=NBII)
summary(nb6a.ex)

eggdata$Male4 <- eggdata$Male.age!=4
nb6d <- gamlss(Eggs~Copulation.duration.seconds + factor(Male4) + Copulation.duration.seconds:factor(Male4) + factor(Female.age), data=eggdata, family=NBII)
LR.test(nb6d, nb6a)

