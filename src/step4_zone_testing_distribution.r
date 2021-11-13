library(circular)

## data
nameExp = "NPA_stich_3"
str_zone = "zoneA"

## parameters
σ2_flock = 2
σ2_stream = 1

## loading
folderDataZone = "../data_tracking/step3_data_zone/"
nameFile = paste0(folderDataZone,nameExp,"/df_velocity/",nameExp,"_",str_zone,".csv")
df = read.csv(nameFile, header = TRUE)
θ = -df$θ

parameter_Wrapped_Gauss <- function(θ,fixedSigma2=-1){
    ## estimating parameter wrapped Gaussian
    ##     f(θ) = C*exp(-(θ-θ_bar)^2/2σ^2)
    ## 1) init
    N = length(θ)
    z = cos(θ) + 1i*sin(θ)
    ## 2) mean/second moment
    z_bar = mean(z)
    R2 = Mod(z_bar)^2
    ## 3) done
    θ_bar = Arg(z_bar)
    ##Re2 = N/(N-1)*(R2-1/N)
    Re2 = R2
    σ2 = log(1/Re2)
    if (fixedSigma2>0){
        σ2 = fixedSigma2
    }
    ## log-likelihood
    #li = dwrappedg(θ,"norm",mean=θ_bar,sd=sqrt(σ2),K=100)
    li = dwrappednormal(θ, mu=circular(θ_bar), sd=sqrt(σ2))
    lli = sum(log(li))
    ## return
    return( c(θ_bar,σ2,lli) )
}

## parameter estimation
θ_tp = seq(-pi,pi,.1)
## Flock
tp = parameter_Wrapped_Gauss(θ,σ2_flock)
θ_bar = tp[1];σ2 = tp[2];lli_flock = tp[3]
f_flock = dwrappednormal(θ_tp, mu=circular(tp[1]), sd=sqrt(σ2_flock))
## Stream
tp = parameter_Wrapped_Gauss(2*θ,4*σ2_stream)
θ_bar = tp[1]/2;σ2 = tp[2]/4;lli_stream = tp[3]
f_stream = dwrappednormal(2*θ_tp, mu=circular(2*θ_bar), sd=sqrt(4*σ2_stream))
## Swarm
N = length(θ)
lli_unif = N*log(1/(2*pi))

## log-likelihood
##---------------
cat("\n",nameExp,str_zone,"\n","\n")
cat("--- log likelihood ---","\n",
    " lli flock  = ",lli_flock,"\n",
    " lli stream = ",lli_stream,"\n",
    " lli swarm  = ",lli_unif,"\n")
cat("--- AIC ---","\n",
    " AIC flock  = ",2*1 - 2*lli_flock,"\n",
    " AIC stream = ",2*1 - 2*lli_stream,"\n",
    " AIC swarm  = ",0 - 2*lli_unif,"\n")
cat("--- BIC ---","\n",
    " BIC flock  = ",1*log(N) - 2*lli_flock,"\n",
    " BIC stream = ",1*log(N) - 2*lli_stream,"\n",
    " BIC swarm  = ",0 - 2*lli_unif,"\n")
## Akaike weights
AIC_all = c(2 - 2*lli_flock,
            2 - 2*lli_stream,
            0 - 2*lli_unif)
delta <- AIC_all - min(AIC_all)
L <- exp(-0.5 * delta) 
AW <- L/sum(L)
cat("--- AW ---","\n",
    " AW flock  = ",AW[1],"\n",
    " AW stream = ",AW[2],"\n",
    " AW swarm  = ",AW[3],"\n")

## plot
plot(θ_tp,f_flock,type='l',lwd=2,col="red",
     xlab="theta",ylab="distribution",cex.lab=1.5,
     xlim=c(-pi,pi),ylim=c(0,.5))
grid()
lines(θ_tp,f_stream,type='l',lwd=3,col="gold")
tp = density(c(θ-2*pi,θ,θ+2*pi),bw=.1) # lines(density(θ,bw=.1)) -> Kernel Density Estimation
lines(tp$x,3*tp$y)
lines(θ_tp,θ_tp*0+1/(2*pi),lwd=2,col="blue")
#lines(θ,θ*0,type='p')
lines(θ[seq(1,23000,by=10)],0*θ[seq(1,23000,by=10)],type='p')
legend("topright",c("flock","stream","swarm"),col=c("red","gold","blue"),lwd=2,cex=1.5)
#dev.print(pdf, 'filename.pdf')
