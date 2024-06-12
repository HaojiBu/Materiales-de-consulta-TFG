#install.packages("MASS")
library(MASS)
#install.packages("lme4")
library(lme4)
#install.packages("pbkrtest")
library(pbkrtest)
library(lmerTest)


# resultados del TFG
nreplicas = 1000
resultadostfg <- matrix(data = NA, nrow = nreplicas, ncol = 4)




# simulación de modelos multinivel con dos factores aleatorios cruzados
# Simulación de parámetros

nsujetos <- 20 #Número de Sujetos del experimento
ntareas <- 6 #Número de tareas

#Intersección y Pendiente del modelo, respectivamente
gamma00 <- 0 
gamma10 <- 0

# y = 

# EFECTOS ALEATORIOS DE LOS SUJETOS

ea_mediasuj <- c(0,0)
var_suj <- 0.3
var_pend_suj <- 0.1
r_sujetos <- 0.5 #Esta es la correlación entre efectos aleatorios de intersección y pendiente de los participantes. Es el parámetro que se modificará para incidir en las covarianzas
covsuj <- sqrt(var_suj) * sqrt (var_pend_suj) * r_sujetos # esto viene de una correlación de -0.5 entre las intersecciones y las pendientes

#MATRIZ VARIANZAS COVARIANZAS DE LOS EFECTOS ALEATORIOS DE LOS SUJETOS
SIGMASUJ <- matrix(data=NA,nrow=2,ncol=2)
SIGMASUJ[1,1]=var_suj
SIGMASUJ[1,2]=covsuj
SIGMASUJ[2,1]=covsuj
SIGMASUJ[2,2]=var_pend_suj

# EFECTOS ALEATORIOS DE LAS TAREAS

ea_mediatarea <- c(0,0)
var_tareas <- 0.2
var_pend_tarea <- 0.1
r_tarea <- 0.5 #Esta es la correlación entre efectos aleatorios de intersección y pendiente de los estímulos o ítems. También se modificará para otras condiciones.
covtarea <- sqrt(var_tareas)* sqrt(var_pend_tarea) * r_tarea 

SIGMATAREA <- matrix(data=NA,nrow=2,ncol=2)
SIGMATAREA[1,1]=var_tareas
SIGMATAREA[1,2]=covtarea
SIGMATAREA[2,1]=covtarea
SIGMATAREA[2,2]=var_pend_tarea


for (k in 1:nreplicas){ # empiezan las iteraciones
  
  
  #mvrnorm es la normal multivariada
  efsujetos <- mvrnorm(nsujetos,ea_mediasuj,SIGMASUJ) #simular la distribución de los efectos aleatorios del sujeto
  
  
  #mvrnorm es la normal multivariada
  eftarea <- mvrnorm(ntareas,ea_mediatarea,SIGMATAREA) #simular la distribución los efectos aleatorios del sujeto
  
  
  #varianza residual de nivel 1
  var_e <- 0.3
  
  y <- c()
  ywide <- matrix(data=NA,nrow = nsujetos, ncol = 2*ntareas)
  idsuj <- c()
  idtarea <- c()
  cont <- 0
  COND <- c(rep(0,ntareas),rep(1,ntareas))
  efsuj1 <- c() # para obtener el efecto de las intersecciones del sujeto
  efsuj2 <- c() # para obtener el efecto de las pendientes del sujeto
  efit1 <- c() # para obtener el efecto de las intersecciones del item
  efit2 <- c() # para obtener el efecto de las pendientes del item
  error <- c()
  dvar <- c()
  for (i in 1:nsujetos){
    cte = 0
    for (j in 1:(2*ntareas)){
      cont = cont + 1
      if (j > ntareas){cte = -ntareas}
      erroraux <- rnorm(1, 0, var_e)
      y[cont] <- gamma00 + gamma10*COND[j] + efsujetos[i,1] + eftarea[cte+j,1] + efsujetos[i,2]*COND[j] +
        eftarea[cte+j,2]*COND[j] + erroraux
      idsuj[cont] <- i
      idtarea[cont] <- cte+j
      ywide[i,j] <- y[cont]
      efsuj1[cont] = efsujetos[i,1]
      efsuj2[cont] = efsujetos[i,2]*COND[j]
      efit1[cont] = eftarea[cte+j,1]
      efit2[cont] = eftarea[cte+j,2]*COND[j]
      error[cont] <- erroraux
      dvar[cont] <- gamma10*COND[j]
    }
  }
  
  aux <- rep(COND,nsujetos)
  COND <- aux
  dlong <- data.frame(idsuj,idtarea,COND, y, dvar,  efsuj1, efsuj2, efit1,efit2, error)
  dwide <- data.frame(ywide)
  #View(dwide)
  dwide$pre <- rowMeans(dwide[, 1:ntareas], na.rm = TRUE)
  dwide$post <- rowMeans(dwide[, (ntareas+1):(2*ntareas)], na.rm = TRUE)
  
  ####################################################################################################
  
  
  ###############                  RESULTADOS              ###########################################
  
  # aproximación clásica
  trelacionadas <- t.test(dwide$pre, dwide$post, paired = TRUE)
  resultadostfg[k,1] <- trelacionadas$p.value
  if (resultadostfg[k,1] < 0.05){resultadostfg[k,3]=1}
  if (resultadostfg[k,1] >= 0.05){resultadostfg[k,3]=0}
  
  # aproximación multinivel con dos factores aleatorios cruzados
  multinivel <- lmer(y ~ COND + (COND|idsuj) + (COND|idtarea), data = dlong)
  
  resultadostfg[k,2] <- summary(multinivel)$coef[2,5]
  if (resultadostfg[k,2] < 0.05){resultadostfg[k,4]=1}
  if (resultadostfg[k,2] >= 0.05){resultadostfg[k,4]=0}
  
  
} # cierra k = nreplicas

#Tasa de rechazos por el modelo clásico y por el modelo multinivel
 
mean(resultadostfg[, 3])
mean(resultadostfg[, 4])
