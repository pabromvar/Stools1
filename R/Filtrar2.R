#' Rule out false positives from data generated with Scalpel and wANNOVAR
#'
#' @param Archivo1
#' @param Archivo2
#' @param NReads
#' @param CovRatio
#' @return Results and False positives excel files 
#' @examples
#' ruleout2()
#' @export

ruleout2= function(Archivo1,Archivo2,CovRatio,NReads){

  Muestra1 = read.csv(Archivo1,header=TRUE, sep = "\t")
  colnames(Muestra1)[colnames(Muestra1)=="Otherinfo.7"] <- "Otherinfo.71"
  Muestra2 = read.csv(Archivo2,header=TRUE, sep = "\t")
  colnames(Muestra2)[colnames(Muestra2)=="Otherinfo.7"] <- "Otherinfo.72"
  
  
  Muestra.1 = transform(Muestra1, Otherinfo.71 = colsplit(Otherinfo.71, split = "\\;", names = c('a', 'b','c','d','e','f', 'g','h','i')))
  a=Muestra.1$Otherinfo.71
  colnames(a)[colnames(a)=="e"] <- "e1"
  Muestra.2 = transform(Muestra2, Otherinfo.72 = colsplit(Otherinfo.72, split = "\\;", names = c('a', 'b','c','d','e','f', 'g','h','i')))
  b=Muestra.2$Otherinfo.72
  colnames(b)[colnames(b)=="e"] <- "e2"
  
  attach(Muestra.1)
  attach(Muestra.2)
  
  
  Datos1 = transform(a, e1 = colsplit(e1, split = "\\=", names = c('Nombre','Cov.Ratio')))
  Datos2 = transform(b, e2 = colsplit(e2, split = "\\=", names = c('Nombre','Cov.Ratio')))
  
  attach(Datos1)
  attach(Datos2)
  
  
  Cov.Ratios1 = e1["Cov.Ratio"]
  Cov.Ratios2 = e2["Cov.Ratio"]
  
  
  genes1 = Muestra1[7]
  genes2 = Muestra2[7]
  
  Comienzo1 = Muestra1[2]
  Comienzo2 = Muestra2[2]
  
  Fin1 = Muestra1[3]
  Fin2 = Muestra2[3]  
  
  Cromosoma1 = Muestra1[1]
  Cromosoma2 = Muestra2[1]
  
  
  Reads.1= transform(Muestra1, Reads1 = colsplit(Otherinfo.9, split = "\\:", names = c('Heterocigosidad', 'Con Indel','Totales')))
  Reads.2= transform(Muestra2, Reads2 = colsplit(Otherinfo.9, split = "\\:", names = c('Heterocigosidad', 'Con Indel','Totales')))
  
  attach(Reads.1)
  attach(Reads.2)
  
  Mutacion1 = Muestra1[4:5]
  Mutacion2 = Muestra2[4:5]
  
  
  ExonicFunc1 = Muestra1[9]
  ExonicFunc2 = Muestra2[9]
  
  Frecuencias1 = Muestra1[11:24]
  Frecuencias2 = Muestra2[11:24]
  
  dbSNP1 = Muestra1[30]
  dbSNP2 = Muestra2[30]
  
  
  tabla1=cbind(Cromosoma1, Comienzo1, Fin1,Mutacion1, ExonicFunc1,genes1,Frecuencias1,dbSNP1,Reads1.Totales,Cov.Ratios1,Reads1.Heterocigosidad)
  tabla2=cbind(Cromosoma2, Comienzo2, Fin2,Mutacion2, ExonicFunc2,genes2,Frecuencias2,dbSNP2,Reads2.Totales,Cov.Ratios2,Reads2.Heterocigosidad)
  
  
  
  masCR.1 = subset(tabla1, Cov.Ratios1 >= CovRatio & Reads1.Totales >= NReads)
  masCR.2 = subset(tabla2, Cov.Ratios2 >= CovRatio & Reads2.Totales >= NReads)
  
  menosCR.1 = subset(tabla1, Cov.Ratios1 < CovRatio | Reads1.Totales < NReads)
  menosCR.2 = subset(tabla2, Cov.Ratios2 < CovRatio | Reads2.Totales < NReads)

  wb = createWorkbook(type="xlsx")
  sheet1 = createSheet(wb, sheetName = Archivo1)
  sheet2 = createSheet(wb, sheetName = Archivo2)
  
  addDataFrame(masCR.1,sheet1)
  addDataFrame(masCR.2,sheet2)
  
  saveWorkbook(wb, "Results.xlsx")
  
  wb2 = createWorkbook(type="xlsx")
  sheet1.1 = createSheet(wb2, sheetName = Archivo1)
  sheet2.2 = createSheet(wb2, sheetName = Archivo2)
  
  addDataFrame(menosCR.1,sheet1.1)
  addDataFrame(menosCR.2,sheet2.2)
  saveWorkbook(wb2, "False Positives.xlsx")
  
  
  devolver = "Done"
  return(devolver)
}
