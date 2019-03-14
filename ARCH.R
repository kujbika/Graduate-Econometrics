library( rugarch ) # Ez az egyváltozós ARCH/GARCH modellkört kezeli, a többváltozósakhoz a rmgarch csomag a jó
#ford closing price based on Eric Zivot: Modeling financial time series
ford <- quantmod::getSymbols( "F", src = "yahoo", from = "1984-02-01", to = "1991-12-31", auto.assign = FALSE )
plot( ford$F.Close )
ford.s <- quantmod::dailyReturn( ford$F.Close ) #logreturn
plot( ford.s )
acf( ford.s ) #no significant autocorrelation among logreturns -> no need for arma model
acf( ford.s^2 ) #we found autocorrelation among the second momentum valued -> need for modeling
acf( abs( ford.s ) ) #same effect

#this means that the vol is clusterized - conditional heteroskedasticity model

hp <- quantmod::getSymbols( "HPQ", src = "yahoo", from = "1984-02-01", to = "1991-12-31", auto.assign = FALSE )
hp.s <- quantmod::dailyReturn( hp$HPQ.Close )

# simulate an ARCH process in first order
SimRes <- ugarchpath( ugarchspec( variance.model = list( model = "sGARCH", garchOrder = c( 1, 0 ) ),
                        mean.model = list( armaOrder = c( 0, 0 ) ),
                        fixed.pars = list( mu = 0, omega = 0.01, alpha1 = 0.8 ) ), n.sim = 250 )
plot( SimRes )
#1: conditional variance
#2: time series plot
par( mfrow = c( 2, 1 ) )
plot( SimRes, which = 2 )
plot( SimRes, which = 1 )
dev.off()
acf( SimRes@path$residSim^2 )
acf( SimRes@path$sigmaSim^2 )

# Engle ARCH LM-próbája:
#p=13
#embed creates a matrix from all the lagged variables. see embed(1 : 30, 3)
mat <- embed( as.numeric( residuals( lm( ford.s ~ 1 ) ) )^2, 13 )
arch.lm <- summary( lm( mat[ , 1 ] ~ mat[ , -1 ] ) )
arch.lm #simple linear regression output
arch.lm$r.squared * length( resid( arch.lm ) )
1 - pchisq( arch.lm$r.squared * length( resid( arch.lm ) ), df = 12 ) #this is the p value: significant arch effect

#mean model is a constant
fit <- ugarchfit( ugarchspec( variance.model = list( model = "sGARCH", garchOrder = c( 1, 1 ) ),
                              mean.model = list( armaOrder = c( 0, 0 ) ) ), ford.s )
fit
#alpha, beta are significant
#alpha+beta<1 -> stable output
#the sum is almost 1 -> volatility is clusterized strongly
#arch lm test is included

plot( fit )
#compare 5-11 and 4-10. 5-11 is very impressive - the model could reduce the autocorrelation from the squared residuals, it captures the process well
#9 is the qq plot, which is also captured by the garch!

####this is the end of the first lesson
AutocorTest( residuals( fit, standardize = TRUE ), lag = 15 )
AutocorTest( residuals( fit, standardize = TRUE )^2, lag = 15 )
ArchTest( residuals( fit, standardize = TRUE ), lag = 12 )


SimRes <- ugarchpath( ugarchspec( variance.model = list( model = "sGARCH", garchOrder = c( 1, 1 ) ),
                                  mean.model = list( armaOrder = c( 0, 0 ) ),
                                  fixed.pars = list( mu = 0, omega = 0.01, alpha1 = 0.05, beta1 = 0.9 ) ),
                      n.sim = 1000 )
plot( SimRes, which = 1 )
plot( SimRes, which = 2 )
hist( SimRes@path$seriesSim )
qqnorm( SimRes@path$seriesSim )
qqline( SimRes@path$seriesSim )

data( "DowJones30" )
head( DowJones30 )
DowJones30$X.Y..m..d <- as.character( DowJones30$X.Y..m..d )

## EGARCH

fit <- ugarchfit( ugarchspec( variance.model = list( model = "sGARCH", garchOrder = c( 1, 1 ) ),
                              mean.model = list( armaOrder = c( 0, 0 ) ) ),
                  DowJones30$XOM )
fit
plot( fit )

fit <- ugarchfit( ugarchspec( variance.model = list( model = "eGARCH", garchOrder = c( 1, 1 ) ),
                              mean.model = list( armaOrder = c( 1, 1 ), include.mean = TRUE ) ),
                  diff( DowJones30$XOM ) )

fit <- ugarchfit( ugarchspec(
  variance.model = list( model = "sGARCH",
                         garchOrder = c( 1, 1 ) ),
  mean.model = list( armaOrder = c( 1, 1 ),
                     include.mean = TRUE,
                     external.regressors = as.matrix( DowJones30$IBM ) ) ),
  diff( DowJones30$XOM ) )
plot( fit )  

## TGARCH

fit <- ugarchfit( ugarchspec( variance.model = list( model = "fGARCH", submodel = "TGARCH", garchOrder = c( 1, 1 ) ),
                              mean.model = list( armaOrder = c( 0, 0 ) ) ),
                  hp.s )
fit
plot( fit )