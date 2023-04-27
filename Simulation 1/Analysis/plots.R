plotConds <- function(obj, conds, yl, leg_labs, sd=NULL){
  
  # obj[,performance,by=c("mechanism", "pm", "corrPred", "pr", "algorithm")]
  corr.labs <- c("independent variables", "correlated variables")
  names(corr.labs) <- c("0", "1")
  pr.labs <- c("5% of relevant", "20% of relevant")
  names(pr.labs) <- c("0.05","0.2")
  mech.labs <- c("MCAR", "MAR", "MNAR")
  names(mech.labs) <- c("mcar", "mar", "mnar")
  
  obj %>% group_by(mechanism, pm, pr, corrPred, algorithm) %>%
    ggplot(aes(pm,performance,color=algorithm)) +
    geom_point() + geom_line() + 
    geom_ribbon(aes(ymin = performance - 1.96 * objSD[,performance], 
                    ymax = performance + 1.96 * objSD[,performance], fill = algorithm), alpha = .3) +
    facet_grid(mechanism ~ pr+corrPred, labeller=labeller(mechanism=mech.labs,
                                                 corrPred=corr.labs,
                                                 pr=pr.labs)) +
    ylab(yl) + 
    scale_x_continuous(name="Percentage of Missing Values") +
    scale_color_brewer(labels=leg_labs, palette="Dark2")
}
