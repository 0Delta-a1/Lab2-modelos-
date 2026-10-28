# =====================================================================
# LABORATORIO MEA: Colinealidad, Selección, Regularización y Validez
# Ejecutar desde una sesión limpia de R. Semilla fija: 2026.
# =====================================================================

set.seed(2026)
options(width = 100, digits = 4, scipen = 999)

# 1. CREACIÓN AUTOMÁTICA DE CARPETAS DE PROYECTO ----------------------
for (d in c("registros", "resultados", "figuras")) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    cat("Carpeta creada:", d, "\n")
  }
}

# 2. CARGA Y VERIFICACIÓN DE PAQUETES ----------------------------------
paquetes <- c("MASS", "car", "leaps", "glmnet", "mlbench")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan)) stop("Instala primero: ", paste(faltan, collapse = ", "))

suppressPackageStartupMessages({
  library(MASS); library(car); library(leaps); library(glmnet); library(mlbench)
})

# Funciones auxiliares de cálculo
rmse <- function(y, pred) sqrt(mean((y - pred)^2, na.rm = TRUE))

condicion <- function(X) {
  ev <- eigen(cor(as.matrix(X)), symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, .Machine$double.eps)
  list(valores_propios = ev, indices = sqrt(max(ev) / ev),
       numero_condicion = sqrt(max(ev) / min(ev)))
}

bkw <- function(X, centrar = TRUE, constante = FALSE) {
  Z <- as.matrix(X)
  if (centrar) Z <- scale(Z, center = TRUE, scale = FALSE)
  if (constante) Z <- cbind(Intercepto = 1, Z)
  normas <- sqrt(colSums(Z^2))
  if (any(normas == 0)) stop("BKW: columna de norma cero")
  Z <- sweep(Z, 2, normas, "/")
  s <- svd(Z)
  indices <- max(s$d) / pmax(s$d, .Machine$double.eps)
  phi <- sweep(s$v^2, 2, pmax(s$d^2, .Machine$double.eps), "/")
  prop <- phi / rowSums(phi)
  rownames(prop) <- colnames(Z)
  colnames(prop) <- paste0("Dim", seq_len(ncol(prop)))
  list(indices = indices, proporciones = prop)
}

vif_seguro <- function(modelo) {
  if (length(attr(terms(modelo), "term.labels")) < 2)
    return(data.frame(variable = character(), VIF = numeric()))
  z <- car::vif(modelo)
  if (is.matrix(z)) z <- z[, "GVIF"]^(1 / (2 * z[, "Df"]))
  data.frame(variable = names(z), VIF = as.numeric(z), row.names = NULL)
}

subconjuntos <- function(datos, respuesta, predictores) {
  f <- reformulate(predictores, response = respuesta)
  rs <- regsubsets(f, data = datos, nvmax = length(predictores),
                   method = "exhaustive", really.big = TRUE)
  s <- summary(rs)
  n <- nrow(datos); sse <- s$rss
  p <- rowSums(s$which)  # Número de parámetros por modelo (incluye intercepto)
  full <- lm(f, data = datos)
  sigma2 <- sum(resid(full)^2) / df.residual(full)
  tab <- data.frame(
    numero_predictores = seq_along(sse), SSE = sse,
    R2_ajustado = s$adjr2,
    Cp = sse / sigma2 - n + 2 * p,
    AIC = n * log(sse / n) + 2 * p,
    BIC = n * log(sse / n) + log(n) * p
  )
  list(ajuste = rs, resumen = s, tabla = tab)
}


# 3. DATOS Y COMPROBACIONES --------------------------------------------
H <- MASS::cement
names(H) <- c("Y", "X1", "X2", "X3", "X4")

vars_cars <- c("MPG.city", "Min.Price", "Price", "Max.Price", "EngineSize",
               "Horsepower", "RPM", "Rev.per.mile", "Fuel.tank.capacity", "Passengers",
               "Length", "Wheelbase", "Width", "Turn.circle", "Rear.seat.room",
               "Luggage.room", "Weight")
C0 <- MASS::Cars93[, vars_cars]
cat("Hald:", dim(H), "\n")
cat("Cars93 original:", dim(MASS::Cars93), "\n")
cat("Filas Cars93 con NA en variables seleccionadas:",
    sum(!complete.cases(C0)), "de", nrow(C0), "\n")
C <- C0[complete.cases(C0), , drop = FALSE]

data("BostonHousing", package = "mlbench")
B <- BostonHousing[, c("medv", "indus", "rm", "ptratio", "lstat")]
B <- B[complete.cases(B), , drop = FALSE]
names(B) <- c("Y", "indus", "rm", "ptratio", "lstat")
stopifnot(nrow(H) == 13, ncol(H) == 5, all(vapply(C, is.numeric, logical(1))))


# 4. F2: DIAGNÓSTICO DE COLINEALIDAD ------------------------------------
cat("\n================ F2: COLINEALIDAD ================\n")
XH <- H[, -1, drop = FALSE]
XC <- C[, setdiff(names(C), "MPG.city"), drop = FALSE]
cat("\nCorrelaciones Hald:\n"); print(round(cor(XH), 3))
corC <- cor(XC)
cat("\nCorrelaciones Cars93:\n"); print(round(corC, 3))

idx <- which(abs(corC) >= .80 & upper.tri(corC), arr.ind = TRUE)
if (nrow(idx)) {
  pares <- data.frame(variable_1 = rownames(corC)[idx[, 1]],
                      variable_2 = colnames(corC)[idx[, 2]], correlacion = corC[idx])
  print(pares[order(-abs(pares$correlacion)), ], row.names = FALSE)
}

mH <- lm(Y ~ ., data = H)
mC <- lm(MPG.city ~ ., data = C)
cat("\nVIF Hald:\n"); print(vif_seguro(mH))
cat("\nVIF Cars93:\n"); print(sort(car::vif(mC), decreasing = TRUE))
cat("\nCondición Hald:\n"); print(condicion(XH))
cat("\nCondición Cars93:\n"); print(condicion(XC))

correr_bkw <- function(X, nombre) {
  configs <- list(centrado_sin_constante = c(TRUE, FALSE),
                  sin_centrar_sin_constante = c(FALSE, FALSE),
                  centrado_con_constante = c(TRUE, TRUE),
                  sin_centrar_con_constante = c(FALSE, TRUE))
  out <- list()
  for (nm in names(configs)) {
    z <- bkw(X, centrar = configs[[nm]][1], constante = configs[[nm]][2])
    out[[nm]] <- z
    cat("\nBKW", nombre, "-", nm, "\n")
    print(data.frame(dimension = seq_along(z$indices),
                     indice_condicion = round(z$indices, 3)))
  }
  invisible(out)
}
bkw_H <- correr_bkw(XH, "Hald")
bkw_C <- correr_bkw(XC, "Cars93")


# 5. F3: SELECCIÓN DE VARIABLES ----------------------------------------
cat("\n================ F3: SELECCIÓN ================\n")
sel_H <- subconjuntos(H, "Y", names(XH))
cat("\nHald, mejor subconjunto por tamaño:\n"); print(sel_H$tabla, row.names = FALSE)

null_H <- lm(Y ~ 1, data = H); full_H <- lm(Y ~ ., data = H)
fwd_H <- stepAIC(null_H, scope = list(lower = null_H, upper = full_H), direction = "forward", trace = FALSE)
bwd_H <- stepAIC(full_H, direction = "backward", trace = FALSE)
both_H <- stepAIC(null_H, scope = list(lower = null_H, upper = full_H), direction = "both", trace = FALSE)
vars_exh_H <- names(coef(sel_H$ajuste, id = which.min(sel_H$tabla$BIC)))[-1]
modelos_H <- list(exhaustivo_BIC = lm(reformulate(vars_exh_H, "Y"), data = H),
                  forward = fwd_H, backward = bwd_H, stepwise = both_H)

sel_B <- subconjuntos(B, "Y", names(B)[-1])
cat("\nBostonHousing, mejor subconjunto por tamaño:\n")
print(sel_B$tabla, row.names = FALSE)
vars_B <- names(coef(sel_B$ajuste, id = which.min(sel_B$tabla$BIC)))[-1]
modelo_B_BIC <- lm(reformulate(vars_B, "Y"), data = B)

pred_C <- setdiff(names(C), "MPG.city")
sel_C <- subconjuntos(C, "MPG.city", pred_C)
cat("\nCars93, mejor subconjunto por tamaño:\n"); print(sel_C$tabla, row.names = FALSE)
vars_BIC_C <- names(coef(sel_C$ajuste, id = which.min(sel_C$tabla$BIC)))[-1]
mC_BIC_exh <- lm(reformulate(vars_BIC_C, "MPG.city"), data = C)
mC_full <- lm(MPG.city ~ ., data = C)
mC_AIC <- stepAIC(mC_full, direction = "both", k = 2, trace = FALSE)
mC_BIC <- stepAIC(mC_full, direction = "both", k = log(nrow(C)), trace = FALSE)


# 6. F4: RIDGE, LASSO Y ELASTIC NET EN HALD ---------------------------
cat("\n================ F4: REGULARIZACIÓN ================\n")
yH <- H$Y; nH <- length(yH); folds_H <- seq_len(nH) 

ajustar_cv <- function(X, y, alpha, foldid) {
  cv.glmnet(as.matrix(X), as.numeric(y), alpha = alpha, foldid = foldid,
            nfolds = length(unique(foldid)), standardize = TRUE,
            type.measure = "mse", grouped = FALSE)
}
cv_r <- ajustar_cv(XH, yH, 0, folds_H)
cv_l <- ajustar_cv(XH, yH, 1, folds_H)
cv_e <- ajustar_cv(XH, yH, .5, folds_H)

res_cv <- function(cv, metodo) data.frame(metodo = metodo,
                                          lambda_min = cv$lambda.min, lambda_1se = cv$lambda.1se,
                                          error_CV_min = min(cv$cvm))
tabla_cv_H <- rbind(res_cv(cv_r, "Ridge"), res_cv(cv_l, "Lasso"),
                    res_cv(cv_e, "Elastic Net alpha=0.5"))
print(tabla_cv_H, row.names = FALSE)

# Guardar Gráficos de Validación Cruzada en 'figuras/'
png("figuras/hald_cv.png", width = 1500, height = 500, res = 140)
par(mfrow = c(1, 3)); plot(cv_r, main = "Ridge"); plot(cv_l, main = "Lasso"); plot(cv_e, main = "Elastic Net")
dev.off()

fit_r <- glmnet(XH, yH, alpha = 0, standardize = TRUE)
fit_l <- glmnet(XH, yH, alpha = 1, standardize = TRUE)
fit_e <- glmnet(XH, yH, alpha = .5, standardize = TRUE)

# Guardar Gráficos de Trayectorias en 'figuras/'
png("figuras/hald_trayectorias.png", width = 1500, height = 500, res = 140)
par(mfrow = c(1, 3))
plot(fit_r, xvar = "lambda", main = "Ridge")
plot(fit_l, xvar = "lambda", main = "Lasso")
plot(fit_e, xvar = "lambda", main = "Elastic Net")
dev.off()

coef_tabla <- function(cv, metodo) {
  a <- as.matrix(coef(cv, s = "lambda.min"))
  b <- as.matrix(coef(cv, s = "lambda.1se"))
  data.frame(variable = rownames(a), metodo = metodo,
             lambda_min = a[, 1], lambda_1se = b[, 1], row.names = NULL)
}
tabla_coef_H <- rbind(coef_tabla(cv_r, "Ridge"), coef_tabla(cv_l, "Lasso"),
                      coef_tabla(cv_e, "Elastic Net"))

# Estabilidad por 100 submuestras
set.seed(2026)
B_est <- 100
freq <- matrix(0, ncol(XH), 2, dimnames = list(colnames(XH), c("Lasso", "ElasticNet")))

for (b in seq_len(B_est)) {
  ids <- sample(seq_len(nH), 10)
  fold_b <- sample(rep(1:5, length.out = length(ids)))
  try({
    l <- cv.glmnet(XH[ids, , drop = FALSE], yH[ids], alpha = 1, foldid = fold_b, nfolds = 5, standardize = TRUE, grouped = FALSE)
    e <- cv.glmnet(XH[ids, , drop = FALSE], yH[ids], alpha = .5, foldid = fold_b, nfolds = 5, standardize = TRUE, grouped = FALSE)
    co_l <- coef(l, s = "lambda.1se")
    co_e <- coef(e, s = "lambda.1se")
    if (!is.null(co_l)) freq[, 1] <- freq[, 1] + (as.matrix(co_l)[-1, 1] != 0)
    if (!is.null(co_e)) freq[, 2] <- freq[, 2] + (as.matrix(co_e)[-1, 1] != 0)
  }, silent = TRUE)
}
freq <- data.frame(variable = rownames(freq), freq / B_est, row.names = NULL)
df_H <- rbind(data.frame(metodo = "Ridge", lambda = fit_r$lambda, df = fit_r$df),
              data.frame(metodo = "Lasso", lambda = fit_l$lambda, df = fit_l$df),
              data.frame(metodo = "Elastic Net", lambda = fit_e$lambda, df = fit_e$df))


# 7. VALIDACIÓN PREDICTIVA CARS93 (100 PARTICIONES 80/20) ------------
cat("\nValidación predictiva Cars93: 100 particiones 80/20\n")
set.seed(2026); B_split <- 100
salidas <- vector("list", B_split)

for (b in seq_len(B_split)) {
  ids <- sample(seq_len(nrow(C)), floor(.8 * nrow(C)))
  train <- C[ids, , drop = FALSE]; test <- C[-ids, , drop = FALSE]
  ols <- lm(MPG.city ~ ., data = train)
  bic <- tryCatch(stepAIC(ols, direction = "both", k = log(nrow(train)), trace = FALSE),
                  error = function(e) ols)
  xt <- as.matrix(train[, pred_C, drop = FALSE])
  xv <- as.matrix(test[, pred_C, drop = FALSE]); yt <- train$MPG.city
  folds <- min(5, length(yt)); foldid <- sample(rep(seq_len(folds), length.out = length(yt)))
  
  r <- cv.glmnet(xt, yt, alpha = 0, foldid = foldid, nfolds = folds, standardize = TRUE, grouped = FALSE)
  l <- cv.glmnet(xt, yt, alpha = 1, foldid = foldid, nfolds = folds, standardize = TRUE, grouped = FALSE)
  e <- cv.glmnet(xt, yt, alpha = .5, foldid = foldid, nfolds = folds, standardize = TRUE, grouped = FALSE)
  
  pred <- c(rmse(test$MPG.city, predict(ols, newdata = test)),
            rmse(test$MPG.city, predict(bic, newdata = test)),
            rmse(test$MPG.city, predict(r, newx = xv, s = "lambda.1se")),
            rmse(test$MPG.city, predict(l, newx = xv, s = "lambda.1se")),
            rmse(test$MPG.city, predict(e, newx = xv, s = "lambda.1se")))
  
  salidas[[b]] <- data.frame(particion = b,
                             metodo = c("MCO completo", "MCO stepAIC-BIC", "Ridge", "Lasso", "Elastic Net"),
                             RMSE = as.numeric(pred))
}

salidas <- do.call(rbind, salidas)
resumen_rmse <- do.call(rbind, lapply(split(salidas$RMSE, salidas$metodo), function(z)
  c(media = mean(z), sd = sd(z))))
resumen_rmse <- data.frame(metodo = rownames(resumen_rmse),
                           RMSE_media = resumen_rmse[, "media"], RMSE_sd = resumen_rmse[, "sd"], row.names = NULL)
resumen_rmse <- resumen_rmse[order(resumen_rmse$RMSE_media), ]

# Guardar Gráfico de RMSE en 'figuras/'
png("figuras/cars93_rmse.png", width = 1100, height = 700, res = 130)
boxplot(RMSE ~ metodo, data = salidas, las = 2, ylab = "RMSE en test",
        main = "Cars93: error en 100 particiones 80/20", xlab = "")
dev.off()


# 8. F5: RESIDUOS E INFLUENCIA -----------------------------------------
cat("\n================ F5: VALIDEZ ================\n")
m_diag <- mC_BIC_exh
n <- nobs(m_diag); p <- length(coef(m_diag))
rs <- rstudent(m_diag); hii <- hatvalues(m_diag); cook <- cooks.distance(m_diag)
dff <- dffits(m_diag); dfb <- dfbetas(m_diag)

casos_cook <- which(cook > 4/n)

# Guardar Gráfico de Diagnósticos en 'figuras/'
png("figuras/cars93_diagnosticos.png", width = 1400, height = 900, res = 130)
par(mfrow = c(2, 2))
plot(fitted(m_diag), rs, xlab = "Ajustados", ylab = "Residuos studentizados", main = "Residuos vs ajustados")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2))
qqnorm(rs, main = "QQ residuos"); qqline(rs)
plot(hii, type = "h", main = "Apalancamiento", ylab = "hii", xlab = "Caso")
abline(h = 2*p/n, lty = 2)
plot(cook, type = "h", main = "Distancia de Cook", ylab = "Cook", xlab = "Caso")
abline(h = 4/n, lty = 2)
dev.off()


# 9. CONSOLIDACIÓN Y EXPORTACIÓN FINAL A CSV Y REGISTROS ---------------
write.csv(salidas, "resultados/cars93_rmse_100_particiones.csv", row.names = FALSE)
write.csv(resumen_rmse, "resultados/cars93_resumen_rmse.csv", row.names = FALSE)
write.csv(sel_H$tabla, "resultados/hald_seleccion.csv", row.names = FALSE)
write.csv(sel_B$tabla, "resultados/boston_seleccion.csv", row.names = FALSE)
write.csv(sel_C$tabla, "resultados/cars93_seleccion.csv", row.names = FALSE)
write.csv(tabla_cv_H, "resultados/hald_cv.csv", row.names = FALSE)
write.csv(tabla_coef_H, "resultados/hald_coeficientes.csv", row.names = FALSE)
write.csv(df_H, "resultados/hald_grados_libertad.csv", row.names = FALSE)
write.csv(freq, "resultados/hald_estabilidad.csv", row.names = FALSE)

writeLines(capture.output(sessionInfo()), "registros/session_info.txt")

cat("\n================ RESUMEN COMPARATIVO FINAL ===============-\n")
if (exists("tabla_coef_H")) {
  cat("\n1. Comparativa de Coeficientes (Hald):\n")
  print(tabla_coef_H, row.names = FALSE)
}

if (exists("resumen_rmse")) {
  cat("\n2. Resumen del Error Predictivo (RMSE en Cars93 - 100 particiones):\n")
  print(resumen_rmse, row.names = FALSE)
}

cat("\n¡PROCESO COMPLETADO EXITOSAMENTE!\n")
cat("- Carpetas automáticas verificadas: 'registros/', 'resultados/', 'figuras/'\n")
cat("- Tablas exportadas en: 'resultados/'\n")
cat("- Gráficos exportados en: 'figuras/'\n")
cat("- Entorno registrado en: 'registros/session_info.txt'\n")