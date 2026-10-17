# BACKBLAZE: APRENDIZAJE ESTADISTICO PARA CONFIABILIDAD
# Descargue este R y backblaze_2013_drive_level.csv en la misma carpeta.
# En RStudio: Session > Set Working Directory > To Source File Location.
# Luego ejecute los bloques en orden con Ctrl+Enter, o todo con Source.
# En una terminal, desde esa carpeta: Rscript ANALISIS_BACKBLAZE.R
#
# Datos publicos: Backblaze Drive Stats 2013; cohorte procesada por disco.
# Fuente: https://www.backblaze.com/cloud-storage/resources/hard-drive-test-data
# No se necesitan los registros diarios originales ni otros archivos R.

# 0. CARPETA Y PAQUETES -------------------------------------------------------
# La carpeta de trabajo debe ser la que contiene estos dos archivos.

# Ejecutar una sola vez si faltan paquetes:
# install.packages(c("survival", "glmnet", "rpart", "randomForestSRC", "gbm", "survivalsvm"))

library(survival)
library(glmnet)
library(rpart)
library(randomForestSRC)
library(gbm)
library(survivalsvm)

# 1. LEER LOS DATOS ----------------------------------------------------------
datos <- read.csv("backblaze_2013_drive_level.csv")
head(datos)

# Una fila por disco. event = 1: falla; event = 0: censura.
# time_days comienza en la primera observacion, no en la instalacion.
datos$first_date <- as.Date(datos$first_date)
datos$capacity_group <- factor(datos$capacity_group)
datos$manufacturer <- factor(datos$manufacturer)

# 2. DESCRIPCION Y KAPLAN-MEIER -----------------------------------------------
resumen <- data.frame(
  Discos = nrow(datos),
  Fallas = sum(datos$event),
  Censuras = sum(datos$event == 0),
  Porcentaje_censura = 100 * mean(datos$event == 0)
)
resumen
summary(datos$time_days)
table(datos$capacity_group, datos$event)

km <- survfit(Surv(time_days, event) ~ capacity_group, data = datos)
colores <- c("#3E4874", "#176D76", "#AB4260", "#5975A4")
plot(km, col = colores, lwd = 2, conf.int = FALSE, mark.time = FALSE,
     xlab = "Dias desde el ingreso", ylab = "Supervivencia", ylim = c(.94, 1))
legend("bottomleft", legend = levels(datos$capacity_group),
       col = colores, lwd = 2, bty = "n")

# 3. ENTRENAMIENTO, VALIDACION Y TEST -----------------------------------------
set.seed(20260704)

# Separar fallas y censuras mantiene proporciones similares en los tres grupos.
fallas <- which(datos$event == 1)
censuras <- which(datos$event == 0)

train_f <- sample(fallas, floor(.60 * length(fallas)))
resto_f <- setdiff(fallas, train_f)
valid_f <- sample(resto_f, floor(.50 * length(resto_f)))
test_f <- setdiff(resto_f, valid_f)

train_c <- sample(censuras, floor(.60 * length(censuras)))
resto_c <- setdiff(censuras, train_c)
valid_c <- sample(resto_c, floor(.50 * length(resto_c)))
test_c <- setdiff(resto_c, valid_c)

datos$grupo <- "test"
datos$grupo[c(train_f, train_c)] <- "train"
datos$grupo[c(valid_f, valid_c)] <- "valid"
table(datos$grupo, datos$event)

# Preparar las covariables.
# Para reproducir los slides se conservan las medianas globales del estudio.
# En un analisis nuevo, estas medianas deben estimarse solo con entrenamiento.
datos$smart_9_first[is.na(datos$smart_9_first)] <-
  median(datos$smart_9_first, na.rm = TRUE)
datos$smart_194_first[is.na(datos$smart_194_first)] <-
  median(datos$smart_194_first, na.rm = TRUE)
datos$smart_5_first[is.na(datos$smart_5_first)] <- 0
datos$smart_197_first[is.na(datos$smart_197_first)] <- 0

datos$log_power_hours <- log1p(datos$smart_9_first)
datos$temp_c <- datos$smart_194_first
datos$log_reallocated <- log1p(datos$smart_5_first)
datos$any_reallocated <- factor(as.integer(datos$smart_5_first > 0))
datos$any_pending <- factor(as.integer(datos$smart_197_first > 0))
datos$first_obs_day <- as.numeric(datos$first_date - min(datos$first_date)) + 1

# Algunos SMART historicos pueden ser posteriores al ingreso.
# La aplicacion es un benchmark, no una validacion prospectiva.
train <- datos[datos$grupo == "train", ]
valid <- datos[datos$grupo == "valid", ]
test <- datos[datos$grupo == "test", ]

formula_surv <- Surv(time_days, event) ~
  capacity_group + manufacturer + log_power_hours + temp_c +
  log_reallocated + any_reallocated + any_pending + first_obs_day

# 4. AFT, COX Y COX PENALIZADO ------------------------------------------------
# 4.1 AFT: Weibull y log-normal, sin penalizacion.
aft_weibull <- survreg(formula_surv, data = train, dist = "weibull")
aft_lognormal <- survreg(formula_surv, data = train, dist = "lognormal")
summary(aft_lognormal)

# Ejemplo: mediana predicha para cinco discos. Puede exceder el seguimiento.
predict(aft_lognormal, newdata = test[1:5, ], type = "quantile", p = .5)

# 4.2 Cox: coeficientes y razones de riesgo.
cox <- coxph(formula_surv, data = train, ties = "efron")
summary(cox)

# 4.3 Matrices numericas para glmnet: sin la columna del intercepto.
x <- model.matrix(formula_surv, data = datos)[, -1]
x_train <- x[datos$grupo == "train", ]
x_valid <- x[datos$grupo == "valid", ]
x_test <- x[datos$grupo == "test", ]
y_train <- Surv(train$time_days, train$event)

# LASSO: alpha = 1. Elegir lambda con el C-index de validacion.
lasso <- glmnet(x_train, y_train, family = "cox", alpha = 1,
                nlambda = 80, standardize = TRUE, cox.ties = "breslow")
lp_valid <- predict(lasso, newx = x_valid, type = "link")
c_lasso <- numeric(ncol(lp_valid))

for (j in seq_len(ncol(lp_valid))) {
  c_lasso[j] <- concordance(Surv(valid$time_days, valid$event) ~ lp_valid[, j],
                            reverse = TRUE)$concordance
}
lambda_lasso <- lasso$lambda[which.max(c_lasso)]
lambda_lasso
coef(lasso, s = lambda_lasso)

# Ridge: alpha = 0. La seleccion sigue el mismo procedimiento.
ridge <- glmnet(x_train, y_train, family = "cox", alpha = 0,
                nlambda = 80, standardize = TRUE, cox.ties = "breslow")
lp_valid <- predict(ridge, newx = x_valid, type = "link")
c_ridge <- numeric(ncol(lp_valid))

for (j in seq_len(ncol(lp_valid))) {
  c_ridge[j] <- concordance(Surv(valid$time_days, valid$event) ~ lp_valid[, j],
                            reverse = TRUE)$concordance
}
lambda_ridge <- ridge$lambda[which.max(c_ridge)]
lambda_ridge

# 5. ARBOL Y RANDOM SURVIVAL FOREST -------------------------------------------
# El arbol usa verosimilitud/desviacion, no el criterio log-rank.
arbol <- rpart(formula_surv, data = train, method = "exp",
               control = rpart.control(cp = .001, minbucket = 80, maxdepth = 5))
arbol$cptable
cp <- arbol$cptable[which.min(arbol$cptable[, "xerror"]), "CP"]
arbol <- prune(arbol, cp = cp)
plot(arbol)
text(arbol, cex = .6)

# RSF: configuracion fija de los slides.
options(rf.cores = 2)
bosque <- rfsrc(formula_surv, data = train, ntree = 150, mtry = 3,
                nodesize = 100, nsplit = 5, importance = "permute",
                seed = 20260704)
rsf_train <- predict(bosque, newdata = train)$predicted
rsf_valid <- predict(bosque, newdata = valid)$predicted
rsf_test <- predict(bosque, newdata = test)$predicted
sort(bosque$importance, decreasing = TRUE)

# 6. BOOSTING COX Y SVM ------------------------------------------------------
# 6.1 Boosting: elegir cuantas etapas conservar, usando validacion.
boost <- gbm(formula_surv, data = train, distribution = "coxph",
             n.trees = 500, interaction.depth = 2, shrinkage = .02,
             n.minobsinnode = 100, bag.fraction = .7,
             train.fraction = 1, verbose = FALSE)

etapas <- c(100, 200, 300, 400, 500)
c_boost <- numeric(length(etapas))

for (j in seq_along(etapas)) {
  riesgo <- predict(boost, newdata = valid, n.trees = etapas[j], type = "link")
  c_boost[j] <- concordance(Surv(valid$time_days, valid$event) ~ riesgo,
                            reverse = TRUE)$concordance
}
mejor_m <- etapas[which.max(c_boost)]
data.frame(Arboles = etapas, C_validacion = c_boost)
mejor_m

# 6.2 SVM: misma submuestra para ambos kernels; solo discos de entrenamiento.
# 100 fallas + 300 censuras. No es una comparacion con igual muestra que RSF.
set.seed(20260721)
indices_svm <- c(sample(which(train$event == 1), 100),
                 sample(which(train$event == 0), 300))

columnas_svm <- all.vars(formula_surv)
svm_train <- train[, columnas_svm]
svm_valid <- valid[, columnas_svm]
svm_test <- test[, columnas_svm]
subtrain <- svm_train[indices_svm, ]

svm_aditivo <- survivalsvm(formula_surv, data = subtrain,
                          type = "regression", kernel = "add_kernel",
                          gamma.mu = .1, opt.meth = "quadprog")
svm_lineal <- survivalsvm(formula_surv, data = subtrain,
                         type = "regression", kernel = "lin_kernel",
                         gamma.mu = .1, opt.meth = "quadprog")

# El SVM predice tiempo: cambiar el signo para obtener un score de riesgo.
aditivo_valid <- -as.numeric(predict(svm_aditivo, newdata = svm_valid)$predicted)
lineal_valid <- -as.numeric(predict(svm_lineal, newdata = svm_valid)$predicted)
aditivo_test <- -as.numeric(predict(svm_aditivo, newdata = svm_test)$predicted)
lineal_test <- -as.numeric(predict(svm_lineal, newdata = svm_test)$predicted)

comparacion_svm <- data.frame(
  Kernel = c("Aditivo", "Lineal"),
  Validacion = c(
    concordance(Surv(valid$time_days, valid$event) ~ aditivo_valid, reverse = TRUE)$concordance,
    concordance(Surv(valid$time_days, valid$event) ~ lineal_valid, reverse = TRUE)$concordance),
  Test = c(
    concordance(Surv(test$time_days, test$event) ~ aditivo_test, reverse = TRUE)$concordance,
    concordance(Surv(test$time_days, test$event) ~ lineal_test, reverse = TRUE)$concordance)
)
comparacion_svm

# Seleccionar por validacion, no por test.
modelos_svm <- list(svm_aditivo, svm_lineal)
svm <- modelos_svm[[which.max(comparacion_svm$Validacion)]]

# 7. COMPARAR TODOS LOS MODELOS -----------------------------------------------
# Mayor score = mayor riesgo. En AFT se invierte el predictor de log-tiempo.
# En RSF, predicted es mortalidad acumulada, no una probabilidad.
riesgo_train <- cbind(
  AFT_Weibull = -predict(aft_weibull, newdata = train, type = "lp"),
  AFT_lognormal = -predict(aft_lognormal, newdata = train, type = "lp"),
  Cox = predict(cox, newdata = train, type = "lp"),
  LASSO = as.numeric(predict(lasso, newx = x_train, s = lambda_lasso, type = "link")),
  Ridge = as.numeric(predict(ridge, newx = x_train, s = lambda_ridge, type = "link")),
  Arbol = as.numeric(predict(arbol, newdata = train)),
  RSF = rsf_train,
  Boosting = predict(boost, newdata = train, n.trees = mejor_m, type = "link"),
  SVM = -as.numeric(predict(svm, newdata = svm_train)$predicted)
)
riesgo_valid <- cbind(
  AFT_Weibull = -predict(aft_weibull, newdata = valid, type = "lp"),
  AFT_lognormal = -predict(aft_lognormal, newdata = valid, type = "lp"),
  Cox = predict(cox, newdata = valid, type = "lp"),
  LASSO = as.numeric(predict(lasso, newx = x_valid, s = lambda_lasso, type = "link")),
  Ridge = as.numeric(predict(ridge, newx = x_valid, s = lambda_ridge, type = "link")),
  Arbol = as.numeric(predict(arbol, newdata = valid)),
  RSF = rsf_valid,
  Boosting = predict(boost, newdata = valid, n.trees = mejor_m, type = "link"),
  SVM = -as.numeric(predict(svm, newdata = svm_valid)$predicted)
)
riesgo_test <- cbind(
  AFT_Weibull = -predict(aft_weibull, newdata = test, type = "lp"),
  AFT_lognormal = -predict(aft_lognormal, newdata = test, type = "lp"),
  Cox = predict(cox, newdata = test, type = "lp"),
  LASSO = as.numeric(predict(lasso, newx = x_test, s = lambda_lasso, type = "link")),
  Ridge = as.numeric(predict(ridge, newx = x_test, s = lambda_ridge, type = "link")),
  Arbol = as.numeric(predict(arbol, newdata = test)),
  RSF = rsf_test,
  Boosting = predict(boost, newdata = test, n.trees = mejor_m, type = "link"),
  SVM = -as.numeric(predict(svm, newdata = svm_test)$predicted)
)

# El mismo calculo del C-index para las nueve columnas.
resultados <- data.frame(Modelo = colnames(riesgo_test),
                          Entrenamiento = NA_real_, Validacion = NA_real_, Test = NA_real_)

for (j in seq_len(ncol(riesgo_test))) {
  resultados$Entrenamiento[j] <- concordance(
    Surv(train$time_days, train$event) ~ riesgo_train[, j], reverse = TRUE)$concordance
  resultados$Validacion[j] <- concordance(
    Surv(valid$time_days, valid$event) ~ riesgo_valid[, j], reverse = TRUE)$concordance
  resultados$Test[j] <- concordance(
    Surv(test$time_days, test$event) ~ riesgo_test[, j], reverse = TRUE)$concordance
}
resultados <- resultados[order(-resultados$Test), ]
print(resultados, digits = 3, row.names = FALSE)

# C-index no significa porcentaje de fallas detectadas.
# Entrenamiento es resustitucion; RSF no usa aqui predicciones OOB.
write.csv(resultados, "comparacion_modelos.csv", row.names = FALSE)
write.csv(comparacion_svm, "comparacion_kernels_svm.csv", row.names = FALSE)

# 8. GRUPOS DE RIESGO DEL BOSQUE ----------------------------------------------
test$cuartil <- cut(riesgo_test[, "RSF"],
                    breaks = quantile(riesgo_test[, "RSF"], seq(0, 1, .25)),
                    include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))

km_riesgo <- survfit(Surv(time_days, event) ~ cuartil, data = test)
plot(km_riesgo, col = colores, lwd = 2, conf.int = FALSE, mark.time = FALSE,
     xlab = "Dias desde el ingreso", ylab = "Supervivencia", ylim = c(.90, 1))
legend("bottomleft", legend = levels(test$cuartil), col = colores, lwd = 2, bty = "n")

tabla_riesgo <- data.frame(
  Cuartil = levels(test$cuartil),
  Discos = as.vector(table(test$cuartil)),
  Fallas = as.vector(table(test$cuartil, test$event)[, "1"]),
  Falla_KM_180_dias = 100 * (1 - summary(km_riesgo, times = 180)$surv)
)
tabla_riesgo$Porcentaje_crudo <- 100 * tabla_riesgo$Fallas / tabla_riesgo$Discos
tabla_riesgo
write.csv(tabla_riesgo, "grupos_riesgo_rsf.csv", row.names = FALSE)

# Son cuartiles descriptivos definidos en test, no una politica de reemplazo.
# Los objetos ajustados quedan disponibles: cox, lasso, bosque, boost, svm, etc.

# 9. PREDECIR PARA UN NUEVO DISCO ---------------------------------------------
# Ejemplo hipotetico: Seagate de 3 TB, 20.000 horas de uso, 32 grados,
# 10 sectores reasignados y 1 sector pendiente. No es un disco del entrenamiento.
# El dia de ingreso 1 corresponde al inicio de la cohorte historica de 2013.
# Estas caracteristicas deben conocerse al momento de realizar la prediccion.
nuevo <- data.frame(
  capacity_group = factor("3 TB", levels = levels(train$capacity_group)),
  manufacturer = factor("Seagate", levels = levels(train$manufacturer)),
  log_power_hours = log1p(20000),
  temp_c = 32,
  log_reallocated = log1p(10),
  any_reallocated = factor(1, levels = levels(train$any_reallocated)),
  any_pending = factor(1, levels = levels(train$any_pending)),
  first_obs_day = 1
)
print(nuevo)

# El bosque ya esta ajustado: solamente se aplica al nuevo perfil.
pred_nuevo <- predict(bosque, newdata = nuevo)
tiempos <- pred_nuevo$time.interest
supervivencia_nuevo <- as.numeric(pred_nuevo$survival[1, ])

# S(t | x): probabilidad de seguir operativo mas alla de t dias.
# 1 - S(t | x): probabilidad de fallar antes o en ese horizonte.
# La curva es escalonada: usar el ultimo tiempo no posterior al horizonte.
# Elegir horizontes dentro del seguimiento, que llega hasta 266 dias.
horizontes <- c(90, 180)
posiciones <- findInterval(horizontes, tiempos)
prediccion_nuevo <- data.frame(
  Dias_desde_ingreso = horizontes,
  Supervivencia = supervivencia_nuevo[posiciones],
  Probabilidad_falla = 1 - supervivencia_nuevo[posiciones],
  Porcentaje_falla = 100 * (1 - supervivencia_nuevo[posiciones])
)
print(prediccion_nuevo, digits = 4, row.names = FALSE)

# El score permite ordenar riesgos; NO es una probabilidad.
# El percentil compara el score nuevo con los scores del entrenamiento.
score_nuevo <- as.numeric(pred_nuevo$predicted)
resumen_nuevo <- data.frame(
  Score_RSF = score_nuevo,
  Percentil_score_en_train = 100 * mean(rsf_train <= score_nuevo),
  Mediana_dias = tiempos[which(supervivencia_nuevo <= .5)[1]]
)
print(resumen_nuevo, digits = 4, row.names = FALSE)

# Mediana_dias = NA: la curva no alcanza 0.5 dentro del seguimiento.
# No significa que el disco nunca fallara ni permite dar un dia exacto de falla.
plot(c(0, tiempos), c(1, supervivencia_nuevo), type = "s", lwd = 2,
     col = colores[1], ylim = c(0, 1), xlab = "Dias desde el ingreso",
     ylab = "Supervivencia predicha para el nuevo disco")
abline(h = .5, col = "gray", lty = 3)
abline(v = horizontes, col = colores[2], lty = 2)
points(horizontes, supervivencia_nuevo[posiciones], pch = 19, col = colores[2])

write.csv(prediccion_nuevo, "prediccion_nuevo_disco.csv", row.names = FALSE)
# Ejemplo didactico: buen C-index no garantiza probabilidades calibradas.
# El benchmark historico todavia no valida un uso prospectivo u operativo.

