#' stack_pca
#'
#' converts two environmental stacks to PCA values
#'
#' @param env1 Environmental stack for time 1
#' @param env2 Environmental stack for time 2
#' @param expvar Proportion of variance of PCA used to select PCA axes. Default to 1 (all PCA axes)
#'
#' @return three objects
#'  summary: PCA summary
#'  pca.time1: PCA for environmental stack at time 1
#'  pca.time2: PCA for environmental stack at time 2
#'
#' @export
#'
#' @examples
stack_pca=function (env1, env2,expvar=1)
{
  env.val1=values(env1)
  keepers1=which(complete.cases(env.val1))
  nas1=which(!complete.cases(env.val1))
  env.val2=values(env2)
  keepers2=which(complete.cases(env.val2))
  nas2=which(!complete.cases(env.val2))

  pca1=prcomp(env.val1[keepers1,], retx = TRUE, center = TRUE,scale = TRUE)
  sum.pca=summary(pca1)
  sum.pca=sum.pca$importance
  sum.pca=as.data.frame(t(sum.pca))
  n=min(which(sum.pca$`Cumulative Proportion`>=expvar))

  env.pca1=env1[[1:n]]
  for (i in 1:n) {
    thislayer = env.pca1[[i]]
    thislayer[nas1] = NA
    thislayer[keepers1] = pca1$x[, i]
    env.pca1[[i]] = thislayer
  }
  names(env.pca1) = paste0("PC", 1:n)

  pca2 = predict(pca1,env.val2[keepers2,])
  env.pca2 = env2[[1:n]]
  for (i in 1:n) {
    thislayer = env.pca2[[i]]
    thislayer[nas2] = NA
    thislayer[keepers2] = pca2[, i]
    env.pca2[[i]] = thislayer
  }
  names(env.pca2) = paste0("PC", 1:n)

  env.pca1 = setMinMax(env.pca1)
  env.pca2 = setMinMax(env.pca2)

  output = list(summary = pca1,pca.time1 = env.pca1, pca.time2 = env.pca2)
  return(output)
}
