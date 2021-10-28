# sets of commonly used custom functions
graphical_pvalue <- function(x){
  x[x < 0.001] = "***"
  x[x > 0.001 & x < 0.01] = "**"
  x[x > 0.01 & x < 0.05] = "*"
  x[x >= 0.05] = "ns"
  print(x)
}

t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  invisible(t.col)
}
