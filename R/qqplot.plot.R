#' Illustration of qqplot.plot()
#'
#' This function generate the qqplot by plotting expected p-values under the fitted model through simulations against the observed p-values. 
#' @param qqplotdata The QQ plot object got from qqplotdata.simu() function.
#' @param seq_inx numeric; QQdata will be thinned every seq_inx. 
#' @param qqplot.axis numeric; the x- and y-axis limits is set from 0 to "qqplot.axis" in the QQ plot. By default, it is 10.
#' @keywords 
#' @export
#' @examples qqplot.plot(qqplotdata,seq_inx=1,qqplot.axis=10)

qqplot.plot <- function(qqplotdata,seq_inx=1,qqplot.axis=10){
  
  QQdata = data.frame(qqplotdata$QQdata)
  obs_lambda = qqplotdata$observedlambda
  m.lambda = qqplotdata$meanEXPlambda
  l.lambda = qqplotdata$lowEXPlambda
  h.lambda = qqplotdata$highEXPlambda
  
  inx <- seq(1,nrow(QQdata),seq_inx)
  QQdata <- QQdata[inx,]
  
  plot(QQdata$mean_log_exp_pvalues, QQdata$log_obs_pvalues, type="l",xlab=expression(Expected~~-log[10](italic(P)~value)), xlim=c(0,qqplot.axis),ylim=c(0,qqplot.axis),ylab=expression(Observed~~-log[10](italic(P)~value)))
  polygon(c(QQdata$lower,rev(QQdata$upper)),c(QQdata$log_obs_pvalues,rev(QQdata$log_obs_pvalues)),col = "grey75", border = FALSE)
  abline(a=0,b=1 ,col = "gray50",lty=2)
  
  mylabel =  bquote(italic(lambda[obs])== .(formatC(obs_lambda,format="f", digits = 3)))
  text(x = 7.5, y = 2.5, labels = mylabel)
  
  mylabel =  bquote(italic(lambda[fit])== .(formatC(m.lambda,format="f", digits = 3)))
  text(x = 7.5, y = 1.5, labels = mylabel)
  
}
  