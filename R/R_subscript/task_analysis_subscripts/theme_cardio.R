require(ggplot2)
require(ggthemes)
theme_cardio <- function(){
  theme_light()+
  theme(plot.title = element_text(size=14))
}
