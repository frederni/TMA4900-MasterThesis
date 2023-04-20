data(mtcars)
library(ggplot2)
library(gridExtra)

p1 <- ggplot(mtcars, aes(x=hp, y=vs)) + 
  geom_point(alpha=.5) + ggtitle("Binomial regression") +
  stat_smooth(method="glm", se=F, method.args = list(family=binomial))
p2 <- ggplot(mtcars, aes(x=hp, y=vs)) + 
  geom_point(alpha=.5) +
  stat_smooth(method="lm", se=F) + ggtitle("Linear regression")

grid.arrange(p1, p2, ncol=2)
ggsave("../figures/linear-vs-logistic-example.pdf")