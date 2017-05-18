#' Rule out false positives from four files generated with Scalpel and wANNOVAR
#'
#' @param Archivo1
#' @param Archivo2
#' @param Archivo3
#' @param Archivo4
#' @param NReads
#' @param CovRatio
#' @return Results and False positives excel files 
#' @examples
#' ruleout4()
#' @export
Filtrar4= function(Archivo1,Archivo2,Archivo3, Archivo4,CovRatio,NReads){
  library("reshape")
  library("xlsx")
  ##Leer todos los archivos##
  Muestra1 = read.csv(Archivo1,header=TRUE, sep = "\t")
  colnames(Muestra1)[colnames(Muestra1)=="Otherinfo.7"] <- "Otherinfo.71"
  Muestra2 = read.csv(Archivo2,header=TRUE, sep = "\t")
  colnames(Muestra2)[colnames(Muestra2)=="Otherinfo.7"] <- "Otherinfo.72"
  Muestra3 = read.csv(Archivo3,header=TRUE, sep = "\t")
  colnames(Muestra3)[colnames(Muestra3)=="Otherinfo.7"] <- "Otherinfo.73"
  Muestra4 = read.csv(Archivo4,header=TRUE, sep = "\t")
  colnames(Muestra4)[colnames(Muestra4)=="Otherinfo.7"] <- "Otherinfo.74"
  

  ##Sacar los datos que nos interesan##
  
  
  Muestra.1 = transform(Muestra1, Otherinfo.71 = colsplit(Otherinfo.71, split = "\\;", names = c('a', 'b','c','d','e','f', 'g','h','i')))
  a=Muestra.1$Otherinfo.71
  colnames(a)[colnames(a)=="e"] <- "e1"
  Muestra.2 = transform(Muestra2, Otherinfo.72 = colsplit(Otherinfo.72, split = "\\;", names = c('a', 'b','c','d','e','f', 'g','h','i')))
  b=Muestra.2$Otherinfo.72
  colnames(b)[colnames(b)=="e"] <- "e2"
  Muestra.3 = transform(Muestra3, Otherinfo.73 = colsplit(Otherinfo.73, split = "\\;", names = c('a', 'b','c','d','e','f', 'g','h','i')))
  c=Muestra.3$Otherinfo.73
  colnames(c)[colnames(c)=="e"] <- "e3"
  Muestra.4 = transform(Muestra4, Otherinfo.74 = colsplit(Otherinfo.74, split = "\\;", names = c('a', 'b','c','d','e','f', 'g','h','i')))
  d=Muestra.4$Otherinfo.74
  colnames(d)[colnames(d)=="e"] <- "e4"
  
  attach(Muestra.1)
  attach(Muestra.2)
  attach(Muestra.3)
  attach(Muestra.4)
  
  Datos1 = transform(a, e1 = colsplit(e1, split = "\\=", names = c('Nombre','Cov.Ratio')))
  Datos2 = transform(b, e2 = colsplit(e2, split = "\\=", names = c('Nombre','Cov.Ratio')))
  Datos3 = transform(c, e3 = colsplit(e3, split = "\\=", names = c('Nombre','Cov.Ratio')))
  Datos4 = transform(d, e4 = colsplit(e4, split = "\\=", names = c('Nombre','Cov.Ratio')))
  
  attach(Datos1)
  attach(Datos2)
  attach(Datos3)
  attach(Datos4)
  
  Cov.Ratios1 = e1["Cov.Ratio"]
  Cov.Ratios2 = e2["Cov.Ratio"]
  Cov.Ratios3 = e3["Cov.Ratio"]
  Cov.Ratios4 = e4["Cov.Ratio"]
  
  
  genes1 = Muestra1[7]
  genes2 = Muestra2[7]
  genes3 = Muestra3[7]
  genes4 = Muestra4[7]
  
  Comienzo1 = Muestra1[2]
  Comienzo2 = Muestra2[2]
  Comienzo3 = Muestra3[2]
  Comienzo4 = Muestra4[2]
  
  Fin1 = Muestra1[3]
  Fin2 = Muestra2[3]  
  Fin3 = Muestra3[3]  
  Fin4 = Muestra4[3]  
  
  Cromosoma1 = Muestra1[1]
  Cromosoma2 = Muestra2[1]
  Cromosoma3 = Muestra3[1]
  Cromosoma4 = Muestra4[1]
  
  
  Reads.1= transform(Muestra1, Reads1 = colsplit(Otherinfo.9, split = "\\:", names = c('Heterocigosidad', 'Con Indel','Totales')))
  Reads.2= transform(Muestra2, Reads2 = colsplit(Otherinfo.9, split = "\\:", names = c('Heterocigosidad', 'Con Indel','Totales')))
  Reads.3= transform(Muestra3, Reads3 = colsplit(Otherinfo.9, split = "\\:", names = c('Heterocigosidad', 'Con Indel','Totales')))
  Reads.4= transform(Muestra4, Reads4 = colsplit(Otherinfo.9, split = "\\:", names = c('Heterocigosidad', 'Con Indel','Totales')))
  
  attach(Reads.1)
  attach(Reads.2)
  attach(Reads.3)
  attach(Reads.4)
  
  Mutacion1 = Muestra1[4:5]
  Mutacion2 = Muestra2[4:5]
  Mutacion3 = Muestra3[4:5]
  Mutacion4 = Muestra4[4:5]
  
  
  ExonicFunc1 = Muestra1[9]
  ExonicFunc2 = Muestra2[9]
  ExonicFunc3 = Muestra3[9]
  ExonicFunc4 = Muestra4[9]
  
  Frecuencias1 = Muestra1[11:24]
  Frecuencias2 = Muestra2[11:24]
  Frecuencias3 = Muestra3[11:24]
  Frecuencias4 = Muestra4[11:24]
  
  dbSNP1 = Muestra1[30]
  dbSNP2 = Muestra2[30]
  dbSNP3 = Muestra3[30]
  dbSNP4 = Muestra4[30]
  
  
  tabla1=cbind(Cromosoma1, Comienzo1, Fin1,Mutacion1, ExonicFunc1,genes1,Frecuencias1,dbSNP1,Reads1.Totales,Cov.Ratios1,Reads1.Heterocigosidad)
  tabla2=cbind(Cromosoma2, Comienzo2, Fin2,Mutacion2, ExonicFunc2,genes2,Frecuencias2,dbSNP2,Reads2.Totales,Cov.Ratios2,Reads2.Heterocigosidad)
  tabla3=cbind(Cromosoma3, Comienzo3, Fin3,Mutacion3, ExonicFunc3,genes3,Frecuencias3,dbSNP3,Reads3.Totales,Cov.Ratios3,Reads3.Heterocigosidad)
  tabla4=cbind(Cromosoma4, Comienzo4, Fin4,Mutacion4, ExonicFunc4,genes4,Frecuencias4,dbSNP4,Reads4.Totales,Cov.Ratios4,Reads4.Heterocigosidad)
  
  masCR.1 = subset(tabla1, Cov.Ratios1 >= CovRatio & Reads1.Totales >= NReads)
  masCR.2 = subset(tabla2, Cov.Ratios2 >= CovRatio & Reads2.Totales >= NReads)
  masCR.3 = subset(tabla3, Cov.Ratios3 >= CovRatio & Reads3.Totales >= NReads)
  masCR.4 = subset(tabla4, Cov.Ratios4 >= CovRatio & Reads4.Totales >= NReads)
  
  menosCR.1 = subset(tabla1, Cov.Ratios1 < CovRatio | Reads1.Totales < NReads)
  menosCR.2 = subset(tabla2, Cov.Ratios2 < CovRatio | Reads2.Totales < NReads)
  menosCR.3 = subset(tabla3, Cov.Ratios3 < CovRatio | Reads3.Totales < NReads)
  menosCR.4 = subset(tabla4, Cov.Ratios4 < CovRatio | Reads4.Totales < NReads)
  
  
  wb = createWorkbook(type="xlsx")
  sheet1 = createSheet(wb, sheetName = Archivo1)
  sheet2 = createSheet(wb, sheetName = Archivo2)
  sheet3 = createSheet(wb, sheetName = Archivo3)
  sheet4 = createSheet(wb, sheetName = Archivo4)
  
  addDataFrame(masCR.1,sheet1)
  addDataFrame(masCR.2,sheet2)
  addDataFrame(masCR.3,sheet3)
  addDataFrame(masCR.4,sheet4)
  
  saveWorkbook(wb, "Results.xlsx")
  
  wb2 = createWorkbook(type="xlsx")
  sheet1.1 = createSheet(wb2, sheetName = Archivo1)
  sheet2.2 = createSheet(wb2, sheetName = Archivo2)
  sheet3.3 = createSheet(wb2, sheetName = Archivo3)
  sheet4.4 = createSheet(wb2, sheetName = Archivo4)
  
  addDataFrame(menosCR.1,sheet1.1)
  addDataFrame(menosCR.2,sheet2.2)
  addDataFrame(menosCR.3,sheet3.3)
  addDataFrame(menosCR.4,sheet4.4)
  
  saveWorkbook(wb2, "False Positives.xlsx")
  
  devolver = "Done"
  return(devolver)
}
