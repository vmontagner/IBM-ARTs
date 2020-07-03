          ####FUNÇÕES####
          #######################
          library(truncnorm)
          #Funcao que calcula sucesso de adquirir territorio/atrair femeas
          func.R = function(R, B) 
            
          {
            (R ^ B / sum(R ^ B)) #R é o valor do recurso do macho guardiao
          }
          
          #nossa boa e velha amiga funcao invlogit
          invlogit = function(x)
          {
            1 / (1 + exp(-x))
          }
          
          #Calculo de oportunidade para seleção
          I_s = function(x) 
          {
            (sum((x - mean(x))^2)/length(x)) / mean(x)^2
          }
          
          #calcula a intensidade de competição espermática (média harmônica)
          SCi = function(x)
          {
            1/mean(1/x)
          }
          
          #Calcula o número de cópulas que cada furtivo tenta
          NF = function(Vf, gv, R, Rmean=5)
          {
            round(Vf*(1+gv*((R-Rmean)/Rmean)))
          }
          
          threshold.harem = function(
            
                                femeas = 1000, #numero de femeas no sistema
                                machos = 1000, #numero de machos no sistema
                                b, #fator de correlação de promiscuidade entre machos guardiao e femeas. 
                                #Quanto maior b, mais copulas furtivas com femeas de harem grandes
                                B, # fator exponencial de importancia de R na aquisição de territorio e femeas
                                #(quanto maior B, mais importante é R)
                                territorio, #numero de territorios no sistema
                                gv, #gv é a seleção sexual pré nos furtivos (scramble competition)
                                V,  # Vantagem do guardiao na competição espermática durante a determinação de paternidade
                                Vf, #Quantidade de cópulas quando R=Rmean
                                
                                r.sd,#Desvio padrão do R
                                sw.sd, #Desvio padrão do SW
                                #para nomear os csv das simulaçoes
                                r.mean, #R medio
                                sw.mean, #SW medio
                                
                                n.geracao #Por quantas gerações simularemos o sistema
          ) 
          
          {
            territorio = round(territorio) #arredonda o valor de territorio (hipercubo gera valores decimais)
            
            #########################################
            #########################################
            
            ###################################
            for (a in 1:n.geracao)
            {
              
              R = rtruncnorm(machos, a = 0, mean = r.mean, sd = r.sd) #valor do recurso de cada macho no momento da maturação sexual
              
              if (a == 1)
              {
                sw = rnorm(machos, sw.mean, sw.sd) #switchpoint genético que determina qual estratégia o macho irá assumir
                sw.f = rnorm(femeas, sw.mean, sw.sd) #switchpoint genético das fêmeas
              }
              else
              {
                sw = sw.geracao1[1:machos] #sw dos 1000 machos da proxima geração
                sw.f = sw.geracao1[(machos+1):(machos+femeas)] #sw das 1000 femeas da proxima geração
              }
              
              limiar = R - sw #Objeto que indica quais machos atingiram o limiar
              guard.R = subset(R, limiar >= 0) #valores de R dos machos guard
              guard.sw = subset(sw, limiar >= 0) #valores de sw dos machos guard
              furt.R = subset(R, limiar < 0) #valores de R dos machos furt
              furt.sw = subset(sw, limiar < 0) #valores de sw dos machos furt
              #Calcula probabilidade de um macho guardiao conseguir um territorio
              prob.guard = func.R(R = guard.R, B = B) 
              
              if(length(guard.R) <= 0)
              {
                stop("Só existem machos furtivos! O sistema de acasalamento por defesa de recursos ou fêmeas deixou de existir")
              }
              
              #Se tiver MENOS individuos guardiões do que territorios
              if (length(guard.R) <= territorio) 
              {
                s.guard = 1:length(guard.R) #todo guardião ganha um territorio
              }
              else #Se tive MAIS individuos guardiões do que territorios
              {
                #sorteia os individuos guardiões, baseado no seu valor de R
                s.guard = sample(length(guard.R), territorio, prob = prob.guard) 
              }
              
              #cria objeto com o valor de sw dos guardioes com territorio
              s.guard.sw = guard.sw[s.guard] 
              #probabilidade de cada macho com territorio conseguir atrair uma femeas
              prob.terr = func.R(guard.R[s.guard], B) 
              
              #hora de distribuir as femeas pelos territoriais
              #criando um data.frame de femeas com o id da femea e o territorio
              #em que ela mora
              femeasH = data.frame(id = paste("F", 1:femeas, sep=""), 
                                  terr = sample(x = s.guard, size = femeas, prob = prob.terr, replace = TRUE),
                                  stringsAsFactors = FALSE)
              rownames(femeasH) = femeasH$id
              #head(femeasH)
              #femeasH["F100",]
              #mean(tamanhoHarem) #deveria ser igual ao Hmean
              tamanhoHarem = tabulate(femeasH$terr, nbins = length(guard.R))
              
              femeasH$Ha = tamanhoHarem[femeasH$terr]
                
              if (length(furt.R) > 0) 
              {
              #calcula quantas tentativas de copulas cada furtivo realiza
              nf = NF(Vf = Vf, R = furt.R, gv = gv, Rmean = r.mean)
              nf = ifelse(nf < 0, 0, nf)
               
              if(sum(nf) != 0)
              {
                #criando um data.frame de tentativas dos machos furtivos
                tentativas = data.frame(furtivo = rep(1:length(furt.R), nf), femea = NA)
                
                #agora temos que sortear uma femea para cada tentativa do macho furtivo.
                for(i in 1:length(furt.R))
                {
                         tentativas$femea[tentativas$furtivo == i] = 
                    sample(femeasH$id, nf[i], replace=FALSE)
                }
                
                #tentativas$femea = sample(femeasH$id, nrow(tentativas), replace=TRUE)
                
                #head(tentativas)
                
                #registrando o numero de femeas por harem de cada femea que vai
                #ser abordada por um furtivo no data.frame de tentativas
                #femeasH["F853", "Ha"]
                tentativas$Ha = femeasH[tentativas$femea,"Ha"] #gambiarra forte aqui
                #head(tentativas)
                
                #calculando as probabilidades
                femeasH$p = invlogit(b*(scale(rank(femeasH$Ha))))
                tentativas$p = femeasH[tentativas$femea, "p"]
                #sorteando os sucessos
                tentativas$success = rbinom(n = nrow(tentativas), size = 1,
                                            prob = tentativas$p)
    
                #head(tentativas)
                #somando os sucessos por macho furtivo
                success = aggregate(tentativas$success, list(tentativas$furtivo), sum)
    
                #Data.frame com todos os furtivos e que tiveram copulas bem sucedidas
                #Foi necessario fazer isso pro data.frame "todos.fitness" tem o mesmo numero
                #de linhas em todas as colunas
                suc.novo = data.frame(furtivo = 1:length(furt.R), tentativas = 0)
                suc.novo[success$Group.1, "tentativas"] = success$x
                
                #calculando o numero de copulas com machos furtivos por femea
                #soh femeas que foram abordadas por furtivos
                fc = aggregate(tentativas$success, list(tentativas$femea), sum)
                
                #agora eh hora de adicionar o numero de copulas totais de cada
                #femea no data.frame delas.
                femeasH$c = 0
                
                #ai pra quem copulou com furtivo, a gente substitui o zero
                #pelo numero de furtivos
                femeasH[fc[,1],"c"] = fc[,2]
                
                #e soma 1 em tudo pra levar em conta que todo mundo copulou com o territorial
                femeasH$c = femeasH$c+1
                
                  #com tudo isso, podemos calcular o sci por macho territorial,
                #(ja que tem uma coluna do macho territorial no data.frame das femeas)
                #calculando o sci por macho territorial
                #head(femeasH)
                sci.m = aggregate(femeasH$c, list(femeasH$terr), SCi)
               
                #recuperando os tamanhos de harem na mesma ordem
                Hterr = aggregate(femeasH$Ha, list(femeasH$terr),mean)
                
                #calcula o sucesso médio dos furtivos
               
                
                #e com isso podemos calcular o scic
                if(var(Hterr[,2]) == 0 | var(sci.m[,2]) == 0)
                  scic = NA
                else
                  scic = cor(Hterr[,2], sci.m[,2])
                
                
              }
              
              else
              {
                femeasH$c = 1
                suc.novo = data.frame(furtivo = 1:length(furt.R), tentativas = 0)
              }
    
              
              ###
                
                sw.geracao1 = rep(NA, length(sw)+length(sw.f)) #cria o vetor vazio do gene sw da proxima geracao
                #copulas = cbind(terr, cop.furt) #matriz com todas as copulas realizadas na geracao
                #dim(copulas)
                maes = sample(x = femeas, size = length(sw.geracao1), replace = TRUE)
                pai_id = rep(NA, length(sw.geracao1)) #vetor vazio com a identidade geral dos pais
                pai_id_g = rep(NA, length(sw.geracao1)) #vetor vazio com a identidade dos pais guardioes
                pai_id_f = rep(NA, length(sw.geracao1)) #vetor vazio com a identidade dos pais furtivos
                pai_sw = rep(NA, length(sw.geracao1))
                
                
                #table(maes)
                #head(femeasH)
                #femeasH[maes[2],]
               
                for (i in 1:length(maes))
                {
                  
                  if (femeasH[maes[i],"c"] == 1) #se TRUE, femea só copulou com o territorial
                  {
                    sw.geracao1[i] = (guard.sw[femeasH[maes[i],"terr"]]/2) + (sw.f[maes[i]]/2) + rnorm(1, mean = 0, sd = 1)
                    pai_id[i] = femeasH[maes[i],"terr"]
                    pai_id_g[i] = femeasH[maes[i],"terr"]
                    pai_sw[i] = guard.sw[femeasH[maes[i],"terr"]]
                  }
                  else
                  {
                    #existe competicao espermatica
                    
                    #linhas do df "tentativas" que tentaram & conseguiram copular com essa femea
                    furtis = tentativas[which(tentativas[, "femea" ]  == femeasH[maes[i], "id"] &
                                          tentativas$success == 1), ] 
                    
                    mae.um = c(femeasH[maes[i], "terr"], furtis$furtivo)
                    mae.prob = c(femeasH[maes[i], "terr"], furtis$furtivo)
                    
                    mae.prob[1] = V / (V + length(furtis$furtivo)) #Probabilidade do macho guardiao
                    mae.prob[-1] = 1 / (V + length(furtis$furtivo)) #Probabilidades dos machos furtivos
                    pai_i = sample(x = mae.um,
                                   size = 1,
                                   prob = mae.prob) #Sorteia o pai
                    
                    if (pai_i == mae.um[1])
                      #Se o pai for o guardião
                    {
                      sw.geracao1[i] = (guard.sw[pai_i]/2) + (sw.f[maes[i]]/2) + rnorm(1, mean = 0, sd = 1)
                      pai_id[i] = pai_i
                      pai_id_g[i] = pai_i
                      pai_sw[i] = guard.sw[pai_i]
                                  } else 
                    {
                      pai_id[i] = pai_i
                      pai_id_f[i] = pai_i
                      pai_sw[i] = furt.sw[pai_i]
                      sw.geracao1[i] = (furt.sw[pai_i]/2) + (sw.f[maes[i]]/2) + rnorm(1, mean = 0, sd = 1)
                      
                    }
                  }
                }
                #calcula quantos filhotes cada macho teve (conta os machos com 0 filhotes)
                ###nao ta contando os territoriais sem territorio, tem que trocar o nbins pelo
                ###numero total de machos, ou incluir zeros no final
                
                pai_id_g = pai_id_g[!is.na(pai_id_g)] #retira os NA do vetor pais guardioes
                pai_id_f = pai_id_f[!is.na(pai_id_f)] #retira os NA do vetor pais furtivos
                
                if(length(pai_id_f) > 0)
                {
                  ans = c(
                    pai.is = I_s(tabulate(pai_id, nbins = machos)),
                    pai.g.mean = mean(tabulate(pai_id_g, nbins = length(guard.R))),
                    pai.g.is = I_s(tabulate(pai_id_g, nbins = length(guard.R))),
                    pai.f.mean = mean(tabulate(pai_id_f, nbins = (machos - length(guard.R)))),
                    pai.f.is = I_s(tabulate(pai_id_f, nbins = (machos - length(guard.R))))
                  )
                }
    
                else
                {
                  ans = c(
                    pai.is = I_s(tabulate(pai_id, nbins = machos)),
                    pai.g.mean = mean(tabulate(pai_id_g, nbins = length(guard.R))),
                    pai.g.is = I_s(tabulate(pai_id_g, nbins = length(guard.R))),
                    pai.f.mean = NA,
                    pai.f.is = NA
                  )
                }
                
        
                if(length(pai_id_f) > 0)
                {
                  delta.fitness = mean(tabulate(pai_id_g, nbins = length(guard.R))) - mean(tabulate(pai_id_f, nbins = length(furt.R))) #número médio de filhotes de guardiões menos o n° médio de filhotes de furtivos
                }
                
                bastardos = length(pai_id_f)/ (length(pai_id_g) + length(pai_id_f)) #proporção de filhotes bastardos nessa geração
                freq.pop = length(furt.R)/ (length(guard.R) + length(furt.R))
                bast.padr = bastardos/freq.pop
                
                if(a == 1)
                {
                  if(length(pai_id_f) > 0)
                  {
                    df.final = matrix(nrow = n.geracao, ncol = 12, dimnames = list(c(),c("Sneaker Proportion", "Mean Switchpoint", "Switchpoint Variation", "Relative Sneaker Offspring", "Delta Fitness", "Mean SCI", "SCIC", "General IS", "Mean Guardian", "Guardian IS", "Mean Sneaker", "Sneaker IS")))
                    df.final[a,] = c(length(furt.sw) / machos, mean(sw), sd(sw), bast.padr, delta.fitness, mean(sci.m$x), scic, ans)
                  }
                  else
                  {
                    df.final = matrix(nrow = n.geracao, ncol = 12, dimnames = list(c(),c("Sneaker Proportion", "Mean Switchpoint", "Switchpoint Variation", "Relative Sneaker Offspring", "Delta Fitness", "Mean SCI", "SCIC", "General IS", "Mean Guardian", "Guardian IS", "Mean Sneaker", "Sneaker IS")))
                    df.final[a,] = c(length(furt.sw) / machos, mean(sw), sd(sw), NA, NA, mean(sci.m$x), scic, ans)
                  }
                }
                
                if(length(pai_id_f) > 0)
                {
                  df.final[a,] = c(length(furt.sw) / machos, mean(sw), sd(sw), bast.padr, delta.fitness, mean(sci.m$x), scic, ans)
                  
                }
                else
                {
                  df.final[a,] = c(length(furt.sw) / machos, mean(sw), sd(sw), NA, NA, mean(sci.m$x), scic, ans)
                }
                
              }
              else
                #caso só tenha guardião na população
              {
                sw.geracao1 = rep(NA, length(sw)+length(sw.f)) #cria o vetor vazio do gene sw da proxima geracao
                maes = rep(NA, length(sw.geracao1)) #vetor vazio com a mae de cada individuo da proxima geracao
                pai_id = rep(NA, length(sw.geracao1)) #vetor vazio com a identidade geral dos pais
                pai_id_g = rep(NA, length(sw.geracao1))
                pai_id_f = rep(NA, length(sw.geracao1))
                
                
                pai_sw = rep(NA, length(sw.geracao1))
                
                
                maes = sample(x = femeas, size = length(sw.geracao1), replace = TRUE)
                
                for (i in 1:length(maes))
                {
                    sw.geracao1[i] = (guard.sw[femeasH[maes[i],"terr"]]/2) + (sw.f[maes[i]]/2) + rnorm(1, mean = 0, sd = 1)
                    pai_id[i] = femeasH[maes[i],"terr"]
                    pai_id_g[i] = femeasH[maes[i],"terr"]
                    pai_sw[i] = guard.sw[femeasH[maes[i],"terr"]]
                }
    
                            pai_id_g = pai_id_g[!is.na(pai_id_g)] #retira os NA do vetor pais guardioes
                pai_id_f = pai_id_f[!is.na(pai_id_f)] #retira os NA do vetor pais furtivos
                
                ans = c(
                  pai.is = I_s(tabulate(pai_id, nbins = machos)),
                  pai.g.mean = mean(tabulate(pai_id_g, nbins = length(guard.R))),
                  pai.g.is = I_s(tabulate(pai_id_g, nbins = length(guard.R))),
                  NA,NA
                )
                  
                if(a == 1)
                {
                  df.final = matrix(nrow = n.geracao, ncol = 12, dimnames = list(c(), c("Sneaker Proportion", "Mean Switchpoint", "Switchpoint Variation", "Relative Sneaker Offspring", "Delta Fitness", "Mean SCI", "SCIC", "General IS", "Mean Guardian", "Guardian IS", "Mean Sneaker", "Sneaker IS")))
                  df.final[a,] = c(length(furt.sw) / machos, mean(sw), sd(sw), NA, NA, NA, NA, ans)
    
                }
                df.final[a,] = c(length(furt.sw) / machos, mean(sw), sd(sw), NA, NA, NA, NA, ans)
    
              }
              
              if (a == n.geracao)
              {
                  guard.fitness = data.frame(morfo = "G", 
                                             R = guard.R, 
                                             SW = guard.sw, 
                                             Territorio = tabulate(s.guard, nbins = length(guard.R)), 
                                             N_Femeas = tabulate(femeasH$terr, nbins = length(guard.R)),
                                             filhotes = tabulate(pai_id_g, nbins = length(guard.R)))
    
                  if(length(furt.R) > 0)
                  {
                    #tag.sneaker = rep("S", length(furt.R))
                    if(length(pai_id_f > 0))
                    {
                      filhSN = tabulate(pai_id_f, nbins = length(furt.R))
                    }
                    else
                    {
                      filhSN = 0
                    }
                    
                    furt.fitness = data.frame(morfo = "S",
                                              R = furt.R, 
                                              SW = furt.sw,
                                              Territorio = 0,
                                              N_Femeas = suc.novo$tentativas,
                                              filhotes = filhSN
                                              ) 
                    
                    todos.fitness = rbind(guard.fitness, furt.fitness)
                  }
    
                  else
                  {
                    todos.fitness = guard.fitness
                  }
    
              }
            }
            return(list(df.final, todos.fitness))
          }  
          
