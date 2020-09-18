
### YR and YC initiators
## group YC or YR initiators

## YR initiators
#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
yr.function <- function(x){
  for(i in 1:length(x)){
    s1 <- data.frame(x)
    yr.initiators <- subset(s1, s1$x == 'CA' | s1$x == 'TG' | s1$x == 'TA' | s1$x == 'CG')
    ini.frame <- data.frame(Ini = 'YR', Freq = sum(yr.initiators$Freq), Prop = sum(yr.initiators$prop))
  }
  return(ini.frame)
}
## YC initiators
#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
yc.function <- function(x){
  for(i in 1:length(x)){
    s1 <- data.frame(x)
    y.initiators <- subset(s1, s1$x == 'CC' | s1$x == 'TC')
    ini.frame <- data.frame(Ini = 'YC', Freq = sum(y.initiators$Freq), Prop = sum(y.initiators$prop))
  }
  return(ini.frame)
}

yr.initiators <- lapply(seqs.prop, yr.function)
yc.initators <- lapply(seqs.prop, yc.function)
## bind lists
f.bind.list <- function(x){
  yr.table <- do.call("rbind",mapply(cbind, x, "SampleID"=names(x),
                                     SIMPLIFY=F))
  return(yr.table)
}
