weighted.var <- function(x, w, wxbar) {
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    return( (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - wxbar)^2) )
}
