library( lattice )
library( latticeExtra )
library( grid )
library( gridExtra )

### preparing data

GSPC <- quantmod::getSymbols( "^GSPC", from = "1900-01-01", auto.assign = FALSE )
LRM <- data.frame( GSPC$GSPC.Close, abs( quantmod::dailyReturn( GSPC$GSPC.Close ) ) )[ -1, ]
names( LRM )[ 2 ] <- "GSPCabsreturn"
LRM$WN <- rnorm( nrow( LRM ) )
LRM$AR1 <- as.numeric( arima.sim( list( ar = 0.7 ), n = nrow( LRM ) ) )
LRM$RW <- cumsum( rnorm( nrow( LRM ) ) )
LRM$index <- 1:nrow( LRM )
str( LRM )

TSs <- c( "GSPCabsreturn", "WN", "AR1", "RW" )
TSvars <- c( GSPCabsreturn = var( LRM$GSPCabsreturn ), WN = 1, AR1 = 1/( 1-0.7^2 ), RW = Inf )

res <- lapply( TSs, function( ts ) xyplot( as.formula( paste0( ts, "~ index" ) ), data = LRM, type = "l", main = ts ) )
do.call( grid.arrange, res )

#### PROPERTIES OF AUTOCORRELATION
##1:acf
res <- lapply( TSs, function( ts ) xyplot( acf ~ lag, data = acf( LRM[[ ts ]], 200, plot = FALSE ), type = "h", main = ts,
                                           ylim = c( -0.1, 1.1 ) ) )
do.call( grid.arrange, res )
#no surprises in 3 out of 4 acf, what about plot(1,1)? Are those ac-s summable?

##2:logarithmic x,y-scale
res <- lapply( TSs, function( ts ) xyplot( abs( acf ) ~ lag, data = acf( LRM[[ ts ]], 200, plot = FALSE ), type = "l", main = ts,
                                           scales = list( log = 10 )  ) )
do.call( grid.arrange, res )
#the plot(1,1) (i.e absolute return of sp500) is linear on log-log scale -> hyperbolic acf!
#its k^-alfa


##3: logarithmic yscale 
res <- lapply( TSs, function( ts ) xyplot( abs( acf ) ~ lag, data = acf( LRM[[ ts ]], 200, plot = FALSE ), type = "l", main = ts,
                                           scales = list(y = list(log = 10 ) ) ) )
do.call( grid.arrange, res )

##4: subset for asymptotics - check abs return of sp500
res <- lapply( TSs, function( ts ) xyplot( abs( acf ) ~ lag, data = acf( LRM[[ ts ]], 200, plot = FALSE ), subset = lag>10,
                                           type = "l", main = ts, scales = list( log = 10 )  ) )
do.call( grid.arrange, res )
xyplot( acf ~ lag, data = acf( LRM$GSPCabsreturn, 200, plot = FALSE ), subset = lag>10, type = c( "l", "r" ),
        scales = list( log = 10 ) )
reg <- lm( log( acf ) ~ log( lag ), data = data.frame( acf = acf( LRM$GSPCabsreturn, lag.max = 200, plot = FALSE )$acf,
                                                       lag = 0:200 ), subset = lag>10 )
summary(reg)
alpha <- -coef( reg )[ 2 ]
alpha
1-alpha/2
#the acf is summable iff alfa is bigger than one(i.e the slope of ) 
#here alfa equals ~0.416, so the abs acf of sp500 return is not summable



fractal::hurstACVF( LRM[[ "GSPCabsreturn" ]], 1000000 )
###ACF MODELLED WITH AR 
xyplot( acf ~ lag, groups = p,
        data = rbind( data.frame( p = "S&P500", acf = acf( LRM$GSPCabsreturn, lag.max = 200, plot = FALSE )$acf, lag = 0:200 ),
                      do.call( rbind, lapply( c( seq( 10, 50, 10 ), 100 ), function( x )
                        data.frame( p = paste0( "AR, p=", x ),
                                    acf = ARMAacf( ar = ar( LRM$GSPCabsreturn, aic = FALSE, order.max = x )$ar, lag.max = 200 ),
                                    lag = 0:200 ) ) ) ),
        type = c( "h", rep( "l", 6 ) ), auto.key = list( columns = 4, points = FALSE, lines = TRUE ), distribute.type = TRUE )

ar( LRM$GSPCabsreturn ) #what is this??
#AR follows the trend up to order p, but after it becomes very bad - CANT MODEL WITH AR
#AR with order 100 is still not good! 

#### PROPERTIES FOR SPECTRE
#asymptotics: f(w)->C * w^(alfa-1) as w->0. Its finite if alfa >1 i.e acf is summable
res <- lapply( TSs, function( ts ) xyplot( spec ~ freq, data = spectrum( LRM[[ ts ]], plot = FALSE ), type =c("p","r"),
                                           main = ts, scales = list( log = 10 ), subset = freq < 0.1,
                                           panel = function( x, y, ... ) {
                                             panel.xyplot( x, y, ... )
                                             grid.text( paste0( "Slope: ", round( lm( y~x )$coefficients[2], 2 ) ),
                                                        unit( 0.3, "npc"),
                                                        unit( 0.95, "npc" ),
                                                        gp = gpar( col =  trellis.par.get()$superpose.line$col[ 1 ] ) )
                                           } ) )
do.call( grid.arrange, res )
#slope is -0.77, so abs acf is not summable! alfa is almost the same as before (in theory
#they have to equal)
reg <- lm( log( spec ) ~ log( freq ), data = spectrum( LRM$GSPCabsreturn, plot = FALSE )[ c( "freq", "spec" ) ],
           subset = freq < 0.1 )
reg
#least absolute deviation regression (l1 norm based linear regression)
reg2 = L1pack::lad( log( spec ) ~ log( freq ), data = spectrum( LRM$GSPCabsreturn, plot = FALSE )[ c( "freq", "spec" ) ],
             subset = freq < 0.1 )
alpha_from_spectre = reg2$coefficients[2] + 1
alpha
alpha_from_spectre
1-alpha/2 #this is the parameter for arfima

 fractal::hurstSpec( LRM[[ "GSPCabsreturn" ]] )

#### PROPERTY FOR VARIANCE OF MEAN

sm2 <- function( ts, min.nm = 10 ) {
  m <- 1:( floor( length( ts )/min.nm ) )
  data.frame( m = m, sm2 = sapply( m, function( m )
    sum( ( apply( matrix( ts[ 1:( floor( length( ts ) / m ) * m ) ], nr = m ), 2, mean )-
             mean( ts ) )^2 )/( floor( length( ts )/m )-1 ) ) )
}
#variance of mean is s^2/n, if X1,X2,..,XN are iid
#the convergence is a little different in terms of AR1. it is s^2 / (n + constant)
#long memory issue: the variance of mean is not linearly decaying

res <- lapply( TSs, function( ts )
  xyplot( sm2 ~ m, data = sm2( LRM[[ ts ]] ), scales = list( log = 10 ), type = c( "p", "r" ), xlab = "m",
          ylab = parse( text = "s[m]^2" ), main = ts,
          abline = c( as.numeric( log10( TSvars[ ts ] ) ), -1 ), panel = function( x, y, ... ) {
            panel.xyplot( x, y, ... )
            grid.text( paste0( "Slope: ", round( lm( y~x )$coefficients[2], 2 ) ), unit( 0.3, "npc"), unit( 0.95, "npc" ),
                       gp = gpar( col = trellis.par.get()$superpose.line$col[ 1 ] ) )
          } ) )
do.call( grid.arrange, res )

xyplot( sm2 ~ m, data = sm2( LRM[[ "GSPCabsreturn" ]] ), scales = list( log = 10 ), type = "p", xlab = "m", ylab = "sm^2",
        main = "GSPCabsreturn", abline = c( as.numeric( log10( TSvars[ "GSPCabsreturn" ] ) ), -1 ) )

reg <- lm( log( sm2 )~log( m ), data = sm2( LRM[[ "GSPCabsreturn" ]] ) )
reg

alpha <- -coef( reg )[ 2 ]
alpha
1-alpha/2

fractal::hurstBlock( LRM$GSPCabsreturn, "aggvar" )

### R/S PROPERTY (according to Hurst 1951)

rs <- function( ts, min.n.block = 2, min.block.size = 8 ) {
  rscalc <- function( ts ) ( max( cumsum( ts-mean( ts ) ) )-
                               min( cumsum( ts-mean( ts ) ) ) )/sd( ts )*( length( ts )-1 )/length( ts )
  m <- min.block.size:( floor( length( ts )/min.n.block ) )
  do.call( rbind, lapply( m, function( m )
    data.frame( m = m, rs = apply( matrix( ts[ 1:( floor( length( ts )/m )*m ) ], nr = m ), 2, rscalc ) ) ) )
}
#in a log-log scale, its a linear line
res <- lapply( TSs, function( ts )
  xyplot( rs ~ m, data = rs( LRM[[ ts ]] ), scales = list( log = 10 ), type = c( "p", "r" ), xlab = "k",
          ylab = "R/S", main = ts, subset = m>100,
          abline = c( as.numeric( log10( TSvars[ ts ] ) ), -1 ), panel = function( x, y, ... ) {
            panel.xyplot( x, y, ... )
            grid.text( paste0( "Slope: ", round( lm( y~x )$coefficients[2], 2 ) ), unit( 0.3, "npc"), unit( 0.95, "npc" ),
                       gp = gpar( col = trellis.par.get()$superpose.line$col[ 1 ] ) )
          } ) )
do.call( grid.arrange, res )

reg <- lm( log( rs ) ~ log( m ), data = rs( LRM[[ "GSPCabsreturn" ]] ), subset = m>100 )
reg
#alpha, again. They have to equal in theory
alpha <- (1-coef( reg )[ 2 ])*2
alpha
1-alpha/2
#fractal does it more sophisticatedly. 0.67 instead of 0.89
fractal::RoverS( LRM[[ "GSPCabsreturn" ]] )
pracma::hurstexp( LRM[[ "GSPCabsreturn" ]] )

### SELF-SIMILARITY, FRACTAL PROPERTY

ScaleMean <- function( ts, min.n.block = 2000, min.block.size = 2 ) {
  m <- min.block.size:( floor( length( ts )/min.n.block ) )
  do.call( rbind, lapply( m, function( m ) {
    mn <- apply( matrix( ts[ 1:( floor( length( ts )/m )*m ) ], nr = m ), 2, mean )
    data.frame( m = m, y = mn, t = seq_along( mn ) )
  } ) )
}
#self-similarity detected! the number on the top shows the block size where we 
#take the average
xyplot( y~t | as.factor( m ), data = ScaleMean( LRM[[ "GSPCabsreturn" ]] ), scales = list( relation = "free" ), type = "l" )

res <- lapply( 2:8, function( c )
  xyplot( GSPCabsreturn[ seq( 1, to = nrow( LRM ), by = c ) ] ~ index[ seq( 1, to = nrow( LRM ), by = c ) ], data = LRM,
          type = "l", xlab = "t", ylab = "y" ) )
do.call( grid.arrange, res )

### ARFIMA simulating to model long memory time series

#acf functions for fractional integrated times series
res <- lapply( c( -0.4, -0.3, 0, 0.1, 0.3, 0.4 ), function( d )
  xyplot( arfima::tacvfARFIMA( dfrac = d, maxlag = 100 )/arfima::tacvfARFIMA( dfrac = d, maxlag = 100 )[1] ~ 0:100, type = "h",
          xlab = "lag", ylab = "ACF", main = paste0( "d=", d ) ) )
do.call( grid.arrange, res )

res <- lapply( c( -0.4, -0.3, 0, 0.1, 0.3, 0.4 ), function( d )
  xyplot( arfima::arfima.sim( 1000, list( d = d ) ) ~ 1:1000, type = "l", xlab = "t", ylab = "y", main = paste0( "d=", d ) ) )
do.call( grid.arrange, res )

### ARFIMA estimation
temp <- arfima::arfima.sim( 1000, list( dfrac = 0.4 ) )
xyplot( temp ~ 1:1000, type = 'l')
arfima::arfima(temp)
temp <- arfima::arfima.sim(10000, list(phi = c(0.3, 0.2), dint = 1, dfrac = 0.4, theta = 0.5))
xyplot(temp~1:10000, type = 'l') #it looks like this cause its integrated 1.4 times
arfima::arfima(temp, order = c(2,1,1)) #order is p, d_integer, q (ARIMA)

fit <- arfima::arfima(LRM[['GSPCabsreturn']])
summary(fit) #fractional integration is 0.2!!! this means long memory issue is detected
#this yield alfa = 1-d*2
## Model diagnostics (did we model the abs return well?):
acf(resid(fit)$Mode1)
acf(resid(fit)$Mode1^2) 
qqnorm(resid(fit)$Mode1)
qqline(resid(fit)$Mode1)


fit <- arfima::arfima(LRM[['GSPCabsreturn']], order = c(1,0,1))
## Model diagnostics (did we model the abs return well?):
acf(resid(fit)$Mode1)
acf(resid(fit)$Mode1^2) 
qqnorm(resid(fit)$Mode1)
qqline(resid(fit)$Mode1)


fit <- fracdiff::fracdiff( LRM[[ "GSPCabsreturn" ]], nar = 0, nma = 0, drange = c( 0, 0.5 ) )
summary( fit )

fit <- fracdiff::fracdiff( LRM[[ "GSPCabsreturn" ]], nar = 1, nma = 1, drange = c( 0, 0.5 ) )
summary( fit )
