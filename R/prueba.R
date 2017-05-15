#' Esta función coge un df y calcula la media de cada columna
#'
#' @param Archivo
#' @param NR
#' @param CovRatio
#' @return Un nuevo df con la media de cada columna del df original 
#' @examples
#' Verd.Posit()
#' @export

Verd.Posit = function(Archivo,CovRatio,NReads){
  Muestra = read.csv(Archivo,header=TRUE, sep = "\t")
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
  mas025 = subset(tabla, Cov.Ratio >= CovRatio & Reads >= NReads)
  menos025 = subset(tabla, Cov.Ratio < CovRatio | Reads < NReads)
  
  resultado = list(mas025)
  
  write.xlsx(mas025, "masCovRatio.xlsx")
  write.xlsx(menos025, "menosCovRatio.xlsx")
  
  ###Gráfico de Calidad###
  F0= nrow(tabla)
  mas01 = subset(tabla, Cov.Ratio >= 0.15)
  F1 = nrow(mas01)
  mas02 = subset(tabla, Cov.Ratio >= 0.25)
  F2 = nrow(mas02)
  F2
  mas03 = subset(tabla, Cov.Ratio >= 0.35)
  F3 = nrow(mas03)
  mas04 = subset(tabla, Cov.Ratio >= 0.45)
  F4 = nrow(mas04)
  mas05 = subset(tabla, Cov.Ratio >= 0.55)
  F5 = nrow(mas05)
  Filtrado = c(F0,F1,F2,F3,F4,F5)
  Pasos = c("F0","F1","F2","F3","F4","F5")
  grafico = data.frame(cbind(Filtrado,Pasos))
  grafico.1= ggplot(grafico, aes(x= Pasos , y = Filtrado)) + geom_point()
  ggsave(filename="Calidad.tiff", plot=grafico.1)
  
  
  ###Gráfico IndelsVSCromosoma###
  
  Muestra.freq = as.data.frame(table(Muestra[1]))
  cromosomas.x = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
  grafico.freq = ggplot(Muestra.freq, aes(x=Var1, y=Freq, fill=Freq)) + geom_bar(stat="identity") + scale_y_continuous(breaks=seq(0,100,1)) + scale_x_discrete(limits=cromosomas.x) + labs(title="Nº de indels (General) VS Cromosoma", x="Cromosoma", y="Nº de indels")
  ggsave(filename="IndelVSCromosoma.tiff", plot=grafico.freq, width = 15 , height = 10 )
  
  ###Gráfico IndelsVSCromosoma (más 0.25)###
  
  Muestra.freq.1 = as.data.frame(table(mas025[1]))
  cromosomas.x = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
  grafico.freq.1 = ggplot(Muestra.freq.1, aes(x=Var1, y=Freq, fill=Freq)) + geom_bar(stat="identity") + scale_y_continuous(breaks=seq(0,100,1)) + scale_x_discrete(limits=cromosomas.x) + labs(title="Nº de indels (>Cov.Ratio) VS Cromosoma", x="Cromosoma", y="Nº de indels")
  ggsave(filename="IndelVSCromosomaPositivo.tiff", plot=grafico.freq.1, width = 15 , height = 10 )
  
  ###Gráfico IndelsVSCromosoma (menos 0.25)###
  
  Muestra.freq.2 = as.data.frame(table(menos025[1]))
  cromosomas.x = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
  grafico.freq.2 = ggplot(Muestra.freq.2, aes(x=Var1, y=Freq, fill=Freq)) + geom_bar(stat="identity") + scale_y_continuous(breaks=seq(0,100,1)) + scale_x_discrete(limits=cromosomas.x) + labs(title="Nº de indels (<Cov.Ratio) VS Cromosoma", x="Cromosoma", y="Nº de indels")
  ggsave(filename="IndelVSCromosomaNegativo.tiff", plot=grafico.freq.2, width = 15 , height = 10 )
  
  devolver = "Hecho"
  return(devolver)
}