
allSame <- function(x) length(unique(x)) == 1

#data pre-processing
data_preparation <- function(input_data,p,t) {
  
  df <- input_data
  
  df1 <- df[complete.cases(df[,c(1,seq(2,ncol(df),t))]),]
  
  df2 <- df1
  
  for (i in (1:nrow(df1))) {
    for (j in (1:(t-1))) {
      if (! (((sum(is.na(df1[i,seq((2+j),ncol(df1),t)])) == 0) | (sum(is.na(df1[i,seq((2+j),ncol(df1),t)])) == p)))) {
        df2[i,c(seq((2+j),ncol(df1),t))] <- rep(NA,p)
      }
    }
  }
  
  return(df2)
}

#finding the coefficients for the group effect
group_effect <- function(input_data, g, p, t, epsilon) {
  
  df3 <- data_preparation(input_data,p,t)
  
  Z.all <- list()
  H.all <- list()
  H_.all <- list()
  S_.all <- list()
  theta_hat.all <- list()
  mu_hat.all <- list()
  
  j <- 1
  
  for (i in seq(2,ncol(df3),t)) {
    
    df4 <- df3[,c(1,i:(i+(t-1)))]
    ANOVA_RM <- AU.2(df4, eps, Group)
    Z.all[[j]] <- ANOVA_RM$group.params$Z
    H.all[[j]] <- as.matrix(ANOVA_RM$group.params$H)
    S_.all[[j]] <- as.matrix(ANOVA_RM$group.params$R_)
    H_.all[[j]] <- as.matrix(ANOVA_RM$group.params$H_)
    theta_hat.all[[j]] <- ANOVA_RM$group.params$teta
    mu_hat.all[[j]] <- ANOVA_RM$group.params$mu
    j <- j + 1
  }
  
  Z <- do.call(cbind, Z.all)
  
  if ((length(unique(H.all)) == 1) & (length(unique(H_.all)) == 1)) {
    H <- H.all[[1]]
    H_ <- H_.all[[1]]
  } else {print("Ошибка: разные матрицы плана")}
  
  theta_hat <- do.call(cbind, theta_hat.all)
  
  mu_hat <- do.call(cbind, mu_hat.all)
  
  S_L <- eye(nrow(H))
  repeat
  {
    S_L_1 <- S_L
    
    R0 <- data.frame(matrix(NA, nrow = p, ncol = p))
    R1 <- data.frame(matrix(NA, nrow = p, ncol = p))
    
    for (i in 1:p)
    {
      for (j in 1:p)
      {
        R0[i,j] <- t(Z.all[[i]] - H%*%theta_hat.all[[i]]) %*% S_L%*% (Z.all[[j]] - H%*%theta_hat.all[[j]])
        R1[i,j] <- t(Z.all[[i]] - H_%*%mu_hat.all[[i]]) %*%S_L%*% (Z.all[[j]] - H_%*%mu_hat.all[[j]])
      }
    }
    
    R0 <- as.matrix(R0)
    R1 <- as.matrix(R1)
    
    A <- (eigen(solve(R1)%*%R0))$vectors[,which.min(eigen(solve(R1)%*%R0)$values)]
    
    df5 <- data.frame(matrix(NA, nrow = nrow(df3), ncol = 1))
    colnames(df5) <- "group"
    df5$group <- df3$group
    for (i in c(2:(t+1)))
    {
      df5[,i] <- as.matrix(df3[,c(seq(i,ncol(df3),t))])%*%A
    }
    
    ANOVA_RM <- AU.2(df5,eps,Group)
    p.val<- ANOVA_RM$P.values[2]
    
    S_L <- ANOVA_RM$group.params$R_
    
    if (FDist2(S_L,S_L_1) <= epsilon) {break}
  }
  
  R0 <- data.frame(matrix(NA, nrow = p, ncol = p))
  R1 <- data.frame(matrix(NA, nrow = p, ncol = p))
  
  for (i in 1:p)
  {
    for (j in 1:p)
    {
      R0[i,j] <- t(Z.all[[i]] - H%*%theta_hat.all[[i]]) %*% S_L%*% (Z.all[[j]] - H%*%theta_hat.all[[j]])
      R1[i,j] <- t(Z.all[[i]] - H_%*%mu_hat.all[[i]]) %*%S_L%*% (Z.all[[j]] - H_%*%mu_hat.all[[j]])
    }
  }
  
  R0 <- as.matrix(R0)
  R1 <- as.matrix(R1)
  
  A <- (eigen(solve(R1)%*%R0))$vectors[,which.min(eigen(solve(R1)%*%R0)$values)]
  
  df5 <- data.frame(matrix(NA, nrow = nrow(df3), ncol = 1))
  colnames(df5) <- "group"
  df5$group <- df3$group
  for (i in c(2:(t+1)))
  {
    df5[,i] <- as.matrix(df3[,c(seq(i,ncol(df3),t))])%*%A
  }
  
  colnames(df5) <- c("group",paste("Time", 1:t, sep = ""))
  
  ANOVA_RM <- AU.2(df5, eps,Group)
  
  p_value<- ANOVA_RM$P.values[2]
  
  return(list(arm_data=df5, A=A, p_value = p_value))
}

#finding the coefficients for the time and interaction effects
time_and_interaction_effect <- function(input_data, g, p, t, effect) {
  
  df3 <- data_preparation(input_data,p,t)
  
  Y.all <- list()
  H.all <- list()
  Lambda.all <- list()
  
  j <- 1
  
  for (i in seq(2,ncol(df3),t)) {
    df4 <- df3[,c(1,i:(i+(t-1)))]
    ANOVA_RM <- AU.2(df4, eps, Group)
    Y.all[[j]] <- ANOVA_RM$ep$Y
    H.all[[j]] <- as.matrix(ANOVA_RM$ep$plan)
    Lambda.all[[j]] <- as.matrix(ANOVA_RM$ep$Lambda)
    j <- j + 1
  }
  
  Y <- as.matrix(do.call(cbind, Y.all))
  
  if (allSame(H.all) == TRUE) {
    H <- H.all[[1]]
    if (effect == "time") {H_u <- as.matrix(H[,-c(seq(1,ncol(H),g))])}
    else if (effect == "interaction") {H_u <- as.matrix(H[,c(seq(1,ncol(H),g))])}
  } else print("Ошибка: разные матрицы план")
  
  if (allSame(lapply(Lambda.all,function(x) rownames(x) <- c())) == TRUE) {Lambda <- Lambda.all[[1]]} else   print("Ошибка: разные корреляционные матрицы")
  
  theta_hat<- ginv(t(H)%*%ginv(Lambda)%*%H)%*%(t(H)%*%ginv(Lambda)%*%Y)
  
  theta_hat_u <- ginv(t(H_u)%*%ginv(Lambda)%*%H_u)%*%(t(H_u)%*%ginv(Lambda)%*%Y)
  
  R0 <- t(Y - H%*%theta_hat)%*%ginv(Lambda)%*%(Y - H%*%theta_hat)
  R1 <- t(Y - H_u%*%theta_hat_u)%*%ginv(Lambda)%*%(Y - H_u%*%theta_hat_u)
  
  
  A <- sapply((eigen(solve(R1)%*%R0))$vectors[,which.min(eigen(solve(R1)%*%R0)$values)], Re)
  
  df5 <- data.frame(matrix(NA, nrow = nrow(df3), ncol = 1))
  colnames(df5) <- "group"
  df5$group <- df3$group
  
  for (i in c(2:(t+1)))
  {
    df5[,i] <- as.matrix(df3[,c(seq(i,ncol(df3),t))])%*%A
  }
  
  colnames(df5) <- c("group",paste("Time", 1:t, sep = ""))
  
  ANOVA_RM <- AU.2(df5,eps,Group)
  
  if (effect == "time") {p_value <- ANOVA_RM$P.values[3]}
  else if (effect == "interaction") {p_value <- ANOVA_RM$P.values[4]}
  
  return(list(arm_data = df5, A=A, p_value = p_value))
}

#building interaction plot for the found coefficients
interaction_plot <- function(arm_data, grnames, g, t, type) {
  
  df5 <- arm_data
  
  val <- do.call("rbind", lapply(c(1:g), function(y) do.call("rbind", lapply(2:(t+1), function(x) na.omit(df5[df5$group == y,x])))))
  
  ng <- lapply(c(1:g), function(y) sapply(2:(t+1), function(x) length(na.omit(df5[df5$group == y,x]))))
  
  gr <- data.frame(grnames,  stringsAsFactors=FALSE)
  gr$counts <- sapply(ng, function(x) sum(x))
  
  group <- data.frame(do.call(c,sapply(c(1:nrow(gr)), function(x) rep(gr[x,1],gr[x,2]))))
  
  rept <- rbind(rep(c(1:t),g),do.call(c, ng))
  time <- data.frame(do.call(c, sapply(c(1:ncol(rept)), function(x) rep(rept[1,x],rept[2,x]))))
  
  interaction_df <- cbind(val,group,time)
  
  colnames(interaction_df) <- c('value','group','time')
  
  if (type == "plot")
  {
    plot(interaction.plot(interaction_df$time, interaction_df$group, interaction_df$value, xlab = 'Момент времени', ylab = 'Среднее значение <<нового признака>>', trace.label = 'Группа',main = "График взаимодействия"))
  }
  
  else if (type == "ggplot")
  {
    interaction_df$time.cat <- factor(interaction_df$time, labels=paste("t", c(1:t), sep=""))
    df <- with(interaction_df , aggregate(value, list(group=group, time=time.cat), mean))
    df$n <- do.call(c, ng)
    
    lp <- ggplot(data = df, aes(x = time, y = x, group = group, shape = group, colour = group, linetype = group)) + geom_line() + geom_point(size = 3)
    
    plot(lp + scale_colour_discrete(name  ="Группа", breaks=grnames, labels=grnames) +
           scale_shape_discrete(name  ="Группа", breaks=grnames, labels=grnames) +
           scale_linetype_discrete(name  ="Группа", breaks=grnames, labels=grnames) +
           theme_bw() + ylab("Среднее значене <<нового признака>> ") + xlab("Время") +
           ggtitle("График взаимодействия"))
  }
  
  
}

MANOVARM_with_missing_data <- function(input_data, g, p, t, effect = "group", epsilon = 0.00001, plt = TRUE, plot_type = "ggplot", grnames = c("group1","group2")) {
  
  if (effect == "group") {
    results <- group_effect(input_data, g, p, t, epsilon)
  } else 
    if ((effect == "time") | (effect == "interaction")) {
      results <- time_and_interaction_effect(input_data, g, p, t, effect)
    }
  
  if (plt == TRUE) {
    arm_data <- results$arm_data
    interaction_plot(arm_data, grnames, g, t, plot_type)
  }
  
  
  #в качестве результата работы программы, будем выдавать p-value, вектор коэффициентов линейной комбинации и график взаимодействия (при необходимости)
  A <- results$A
  p_value <- results$p_value
  
  return(list(A=A, p_value = p_value))
}
