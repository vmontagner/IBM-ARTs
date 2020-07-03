          parei = 3000  

setwd("~/Área de Trabalho/Simulação_04_07_2020")
            factors <- c("p", "S_G", "T", "g", "V", "C", "R_SD") #nomes dos parametros
                        q = c("qunif")
                        q.arg <- list( 
                          list(min = -3, max=3), #p ator de correlação de promiscuidade entre machos guardiao e femeas
                          list(min= 0, max=5), #S_G fator exponencial de importancia de R na aquisição de territorio e femeas
                          list(min = 100, max = 800), #T = territorio
                          list(min = 0, max = 1), #g
                          list(min = 1/5, max = 5), #V
                          list(min = 0.5, max = 8), #C
                          list(min = 0.25, max = 3) #R_SD
                        )
                        
                        library(pse)
                        myLHS <- LHS(model = NULL, method = "random", factors = factors, N = 3000, q = q, q.arg = q.arg, nboot=50)
            
                        save(myLHS, file = "myLHS.RData")
                        load("myLHS.RData")
                        LHS = get.data(myLHS)
                            geracao = 400 #numero de gerações
                            sw.mean = 5
                            sw.sd = 1
                            r.mean = 5
                            #r.sd = r.mean*0.5
                            
                            #Primeiro, é preciso gerar os valores de LHS no script do HIPERCUBO
                            repeticoes = LHS[rep(1:nrow(LHS), each = 1),] 
                            WID = 1:nrow(repeticoes)
                            filename = paste("Sim_", WID, "_.csv", sep = "")
                            filename3 = paste("Fitness_Sim_", WID, "_.csv", sep="")
                            repeticoes = data.frame(WID, filename, repeticoes)
                            names(repeticoes) = c("WID", "Filename", "p", "S_G", "T", "g", "V", "C", "R_SD")
                            head(repeticoes)
                            write.csv(x = repeticoes, file = "LHS_parametros_entrada.csv")
                            
                            
                            resultados = data.frame(WID)
                            resultados$filename = as.character(filename)
                            resultados$freq.furtivos = NA
                            resultados$limiar.medio = NA
                            resultados$limiar.var = NA
                            resultados$bastardos = NA
                            resultados$Delta.fitness = NA
                            resultados$SCIm = NA
                            resultados$SCIC = NA
                            resultados$pai.is = NA
                            resultados$pai.g.mean = NA
                            resultados$pai.g.is = NA
                            resultados$pai.f.mean = NA
                            resultados$pai.f.is = NA
                            names(resultados) = c("WID", "Filename","Sneaker Proportion", "Mean Switchpoint", "Switchpoint Variation", "Relative Sneaker Offspring", "Delta Fitness", "Mean SCI", "SCIC", "General IS", "Mean Guardian", "Guardian IS", "Mean Sneaker", "Sneaker IS")
                            resultados$Filename = as.character(resultados$Filename)
                            head(resultados)
                            
                            resultados3 = data.frame(WID)
                            resultados3$filename = as.character(filename3)
                            resultados3$morfo = NA
                            resultados3$R = NA
                            resultados3$SW = NA
                            resultados3$Territorio = NA
                            resultados3$N_Femeas = NA
                            resultados3$filhotes = NA
                            names(resultados3) = c("WID", "Filename", "Morph",  "Body Condition", "SW", "Teritory", "N Females", "N offspring")
                            head(resultados3)
                            
                            
                            filename = paste("csv/",filename, sep = "")
                            filename3 = paste("fitness/", filename3, sep = "")
                            
                            for(i in parei:nrow(repeticoes))
                            {
                              if(i == 1)
                              {
                              start.time <- Sys.time()
                              }
                              
                              set.seed(repeticoes[i,1])
                              #roda uma simulacao
                              simi = threshold.harem(
                                b = repeticoes[i,3], 
                                B = repeticoes[i,4], 
                                territorio = repeticoes[i,5], 
                                gv = repeticoes[i,6], 
                                V = repeticoes[i,7], 
                                Vf = repeticoes[i,8], 
                                r.sd = repeticoes[i,9], 
                                sw.sd = sw.sd, 
                                sw.mean = sw.mean, 
                                r.mean = r.mean, 
                                n.geracao = geracao)
                         
                              
                              resultados$`Sneaker Proportion`[i] = simi[[1]][geracao,1]
                              resultados$`Mean Switchpoint`[i] = simi[[1]][geracao,2]
                              resultados$`Switchpoint Variation`[i] = simi[[1]][geracao,3]
                              resultados$`Relative Sneaker Offspring`[i] = simi[[1]][geracao,4]
                              resultados$`Delta Fitness`[i] = simi[[1]][geracao,5]
                              resultados$`Mean SCI`[i] = simi[[1]][geracao,6]
                              resultados$SCIC[i] = simi[[1]][geracao,7]
                        
                              resultados$`General IS`[i] = simi[[1]][geracao,8]
                              resultados$`Mean Guardian`[i] = simi[[1]][geracao,9]
                              resultados$`Guardian IS`[i] = simi[[1]][geracao,10]
                              resultados$`Mean Sneaker`[i] = simi[[1]][geracao,11]
                              resultados$`Sneaker IS`[i] = simi[[1]][geracao,12]
                              
                              resultados3$Morph = simi[[2]][,1]
                              resultados3$`Body Condition` = simi[[2]][,2]
                              resultados3$SW = simi[[2]][,3]
                              resultados3$Teritory = simi[[2]][,4]
                              resultados3$`N Females` = simi[[2]][,5]
                              resultados3$`N offspring`= simi[[2]][,6]
                              
                              write.csv(x = simi[[1]], file = filename[i])
                              write.csv(x = simi[[2]], file = filename3[i])
                        
                              if (i == nrow(repeticoes))
                              {
                                end.time <- Sys.time()
                                time.taken <- end.time - start.time
                                time.taken
                                
                                write.table(resultados, file = "resultados_ultima_geracao.csv", sep = ",", row.names = FALSE)
            
                              }
                              
                              cat("simulacao", i, "concluida. Hold on.\n")
                            }
                                                                                                                