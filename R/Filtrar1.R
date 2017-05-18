#' Rule out false positives from data generated with Scalpel and wANNOVAR
#'
#' @param Archivo1
#' @param NReads
#' @param CovRatio
#' @return Results and False positives excel files 
#' @examples
#' ruleout1()
#' @export

ruleout1 = function(Archivo1,CovRatio,NReads){
  library("reshape")
  library("ggplot2")
  library("xlsx")
  
  Muestra = read.csv(Archivo1,header=TRUE, sep = "\t")
  Muestra.1 = transform(Muestra, Otherinfo.7 = colsplit(Otherinfo.7, split = "\\;", names = c('a', 'b','c','d','e','f', 'g','h','i')))

  attach(Muestra.1)
  Datos = transform(Otherinfo.7, e = colsplit(e, split = "\\=", names = c('Nombre','Cov.Ratio')))
  attach(Datos)
  Cov.Ratios = e["Cov.Ratio"]
  genes = Muestra[7]
  Comienzo = Muestra[2]
  Fin = Muestra [3]
  Cromosoma = Muestra [1]
  Reads.1= transform(Muestra, Reads = colsplit(Otherinfo.9, split = "\\:", names = c('Heterocigosidad', 'Con Indel','Totales')))
  attach(Reads.1)
  Reads= Reads.Totales
  Mutacion = Muestra[4:5]
  
  
  tabla=cbind(Cromosoma, Comienzo, Fin,Mutacion, ExonicFunc.refgene,genes,ExAC_Freq,dbSNP,Reads,Cov.Ratios)
  masCR = subset(tabla, Cov.Ratio >= CovRatio & Reads >= NReads)
  menosCR = subset(tabla, Cov.Ratio < CovRatio | Reads < NReads)
  

  write.xlsx(masCR, "Results.xlsx")
  write.xlsx(menosCR, "False Positives.xlsx")

  devolver = "Done"
  return(devolver)
}

