data(mtcars)
library(ggplot2)
library(cowplot)
p1 <- ggplot(mtcars, aes(x=hp, y=vs)) + 
  geom_point(alpha=.5) + ggtitle("Binomial regression") +
  stat_smooth(method="glm", se=F, method.args = list(family=binomial)) +
  ylim(c(-0.2,1.01))
p2 <- ggplot(mtcars, aes(x=hp, y=vs)) + 
  geom_point(alpha=.5) +
  stat_smooth(method="lm", se=F) + ggtitle("Linear regression") +
  ylim(c(-0.2,1.01))

plot_grid(p1, p2, ncol = 2, align = "v", axis = "tb")
ggsave("../figures/linear-vs-logistic-example.pdf")
