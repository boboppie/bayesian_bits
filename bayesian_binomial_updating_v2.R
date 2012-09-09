sim_bayes<-function(p=0.5,N=10,y_lim=15,prior_a=1,prior_b=1)
{
  success<-0
  curve(dbeta(x,prior_a,prior_b),xlim=c(0,1),ylim=c(0,y_lim),xlab='p',ylab='Posterior Density',lty=2)
  legend('topright',legend=c('Prior','Updated Posteriors','Final Posterior'),lty=c(2,1,1),col=c('black','black','red'))
  for(i in 1:N)
  {
    if(runif(1,0,1)<=p)
      success<-success+1
    
    curve(dbeta(x,success+prior_a,(i-success)+prior_b),add=TRUE)
    print(paste(success,"successes and ",i-success," failures"))
  }
  curve(dbeta(x,success+prior_a,(i-success)+prior_b),add=TRUE,col='red',lwd=1.5)
}

sim_bayes(p=0.2, N=50,prior_a=10,prior_b=10)