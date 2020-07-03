#######################################
#######################################
#PREPARAÇÃO PARA GERAÇÃO DE GRAFICOS###
########################################
#Ir para a pasta e carregar o csv com resultados finais das simulações
saida = read.csv("resultados_ultima_geracao.csv", sep = ",", row.names = 1)
saida = saida[,-1] #Resultados do LHS
head(saida)
entrada = read.csv("LHS_parametros_entrada.csv", sep = ",", row.names = 1)
head(entrada)
entrada = entrada[,-c(1,2)]
#entrada = read.csv("LHS_parametros_entrada.csv", sep = ",", row.names = 1)
#entrada = entrada[1:5000,3:9]

library(sensitivity)
library(pse)

chosenlines = which(complete.cases(saida))
length(chosenlines)

head(saida)
ab = pcc(X = entrada[chosenlines,], y = saida[chosenlines,1], rank=TRUE)
ab$PRCC$original

nome = c("p", "S_G", "Territory", "S_S", "V", "C", "R_SD")
data2 = data.frame(row.names = nome,
                   PropSneaker = ab$PRCC$original,
                   MeanSwitch = NA,
                   VarSwitch = NA,
                   SneakOffsp = NA,
                   D.Fit = NA,
                   MeanSCI = NA,
                   SCIC = NA,
                   IS = NA,
                   M.Off.Guard = NA,
                   ISGuard = NA,
                   M.Off.Sneak = NA,
                   ISSneak = NA)

for (i in 2:ncol(saida))
{
  ab = pcc(X = entrada[chosenlines,], y = saida[chosenlines,i], rank=TRUE)
  data2[,i] = ab$PRCC$original
}

data2 = round(data2, 3)
write.csv(x = data2, file = "PRCC.csv" )
write.csv2(x = data2, file = "PRCC.csv")  
write.table(data2, file= "PRCC_table.csv", sep = ",")
