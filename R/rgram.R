rgram <- function (object, xlim = NULL, ylim = NULL, xlab = "Count", ylab = NULL, main = NULL, 
breaks = NULL, border = "black", fill = "lightgray", col = "blue", 
         lwd = 2, pch = 19, lty = 1,  axes = TRUE, width = NULL, plot = TRUE, 
         style = c("hanging", "standing", "suspended"), scale = c("sqrt", "raw"), 
         max = NULL, ...) 
{
    UseMethod("rgram")
}
