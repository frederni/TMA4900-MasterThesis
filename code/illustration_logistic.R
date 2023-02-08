data(mtcars)
library(ggplot2)
library(gridExtra)

fit0 = lm(vs ~ hp, data=mtcars)
fit1 = glm(vs ~ hp, data=mtcars, family=binomial)

p1 <- ggplot(mtcars, aes(x=hp, y=vs)) + 
  geom_point(alpha=.5) +
  stat_smooth(method="glm", se=FALSE, method.args = list(family=binomial))
p2 <- ggplot(mtcars, aes(x=hp, y=vs)) + 
  geom_point(alpha=.5) +
  stat_smooth(method="lm", se=FALSE)

grid.arrange(p1, p2, ncol=2)
pdf("../logistic_vs_linear.pdf")
dev.off()

#define new data frame that contains predictor variable
  
#use fitted model to predict values of vs
