#load data 50 obs of 3 different types of flower. 4 cont variables
iris
head(iris,3)

# log transform  - box cox 
log.ir <- log(iris[, 1:4])
ir.species <- iris[, 5]

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
ir.pca <- prcomp(log.ir,
                 center = TRUE,
                 scale. = TRUE)

ir.pca
# prcomp data type
# rotations are the loadings of the coefficients of the linear combo of the 4 variable 
print(ir.pca)

# plot variances against variables
# see which variables explain most of the variability in the data (to retain for further analysis)
#clearly 1 and 2
plot(ir.pca, type = "l")

# summary method
# clearly first 2 variables explain 96% of data variance
summary(ir.pca)

# Predict PCs - use if we have new data and we want to predict their PCA values
predict(ir.pca, 
        newdata=tail(log.ir, 2))


library(ggplot2)
#create data frame
scores = as.data.frame(ir.pca$x)
ir.species


scores = cbind(scores,ir.species)
scores
#plot
ggplot(data = scores, aes(x = PC1, y = PC2, color=ir.species)) +
  geom_point(shape=1)+
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  ggtitle("Iris Data")

install.packages("devtools")
library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups = ir.species, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
