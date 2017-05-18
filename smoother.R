test_function = function(x, type="Sine", n=100){
  switch(type,
         Sine = {return(sin(2*pi*x))},
         W = {return(weierstrass_function(x, n))}
         )
}

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
  #select evenly-spaced knots
  K <- seq(min(x), max(x), length.out = (k+2))[2:(k+1)]
  A <- matrix(, length(x), length(K))
  for (m in 1:length(x)){
    for (n in 1:length(K)){
      row[n] <- (max(0,x[m]-K[n]))^l
    }
    A[m,] <- row
  }
  
  X <- cbind(poly(x, l, raw = T), A)
  w_hat <- coef(lm(y~X))
  pred <- cbind(rep(1, length(sq)), cbind(poly(sq, l, raw = T), A)) %*% w_hat 
}

kernel_function = function(t, kernel = "Epanechnikov"){
  if(abs(t)>1){
    return(0)
  }else{
    switch(kernel,
           Uniform = {return(1/2)},
           Triangular = {return(1-abs(t))},
           Epanechnikov = {return(3/4*(1-t^2))},
           Tricube = {return(70/81*(1-abs(t)^3)^3)},
           Gaussian = {return(1/sqrt(2*pi)*exp(-1/2*t^2))}
           )
  }
}

Nadaraya_Watson = function(x, y, sq, lambda=0.2){
  f <- function(sqs){
    for(i in 1:length(sqs)){
      numerator = 0
      denominator = 0
      for (j in 1:length(x)){
        tj <- abs(x[j]-sq[i])/lambda
        kj <- kernel_function(tj)
        kyj <- kj * y[j]
        numerator = numerator + kyj
        denominator = denominator + kj
      }
      if (denominator == 0){stop("k(x0, xi) sums to zero")}
      else{
        yi_hat = numerator/denominator
        y_hat[i] <- yi_hat
      }
    }
    return(y_hat)
  }
  return(f(sq))
}

