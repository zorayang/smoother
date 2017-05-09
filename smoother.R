test_function = function(x, type="Sine", n=100){
  switch(type,
         Sine = {return(sin(2*pi*x))},
         W = {return(weierstrass_function(x, n))}
         )
}

#nth degree Weierstrass function, periodical with 2pi
#references: https://sites.math.washington.edu/~conroy/general/weierstrass/weier.htm
weierstrass_function = function(x, n){
  y = 0
  for (i in 1:n){
    y <- y + (1/2^i)*sin(4*(2^i)*x)
  }
  return(y)
}

data.gen = function(a = 0, b = 1, type="Sine", len=100, noise=.17, n=100){
  x <- seq(a, b, length.out = len)
  switch(type,
         Sine = {y <- sin(2*pi*x) + runif(len, -noise, noise)},
         W = {y <- weierstrass_function(x, n) + runif(len, -noise, noise)}
         )
  return(data <- cbind(x, y))
}

#least sqaure solution
#http://www.sci.ccny.cuny.edu/~szlam/2013-fall-366/least_squares.pdf
poly_approx = function(x, y, sq, n, method = "QR"){
  A <- cbind(rep(1, length(x)), poly(x, n, raw = T))
  b <- y
  AT <- t(A)
  w_hat <- qr.solve(AT %*% A, AT %*% b) # solve for the normal equation
  # coef(lm(y~x))
  pred <- cbind(rep(1, length(sq)), poly(sq, n, raw = T)) %*% w_hat 
  return(pred)
}

poly_approx = function(x, y, sq, n){
  w_hat <- coef(lm(y~poly(x, n, raw = T)))
  pred <- cbind(rep(1, length(sq)), poly(sq, n, raw = T)) %*% w_hat 
  return(pred)
}

#https://faculty.washington.edu/heagerty/Courses/b571/homework/spline-tutorial.q
natural_spline = function(x, y, sq, k, l){
  #select evenly spaced knots
  K <- seq(min(x), max(x), length.out = (k+2))[2:(k+1)]
  
  #there's probably a better way to do this:
  A <- matrix(, length(x), length(K))
  for (m in 1:length(x)){
    row <- length(K)
    for (n in 1:length(K)){
      row[n] <- (max(0,x[m]-K[n]))^l
    }
    A[m,] <- row
  }
  
  X <- cbind(poly(x, l, raw = T), A)
  w_hat <- coef(lm(y~X))
  pred <- cbind(rep(1, length(sq)), cbind(poly(sq, l, raw = T), A)) %*% w_hat 
}