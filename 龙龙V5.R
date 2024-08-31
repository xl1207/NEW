library(broom)
library(readxl)
library(flextable)
library(officer)

mydata <- read_excel("454加34填充好了.xlsx", sheet = "Sheet1")
mydatacopy <-mydata
# 循环遍历自变量列
models <- list()
for (i in 2:ncol(mydatacopy)) {
  y <- data[, 1]
  y<-unlist(y)
  x<- data[, i]
  x<-unlist(x)
  # 进行逻辑回归
  model <- glm(y ~ x, data = data, family = binomial(link = "logit"))
  models[[i]] <- model
}
models <- models[-1]
#造数据框
data <- data.frame(column1 = c("value1", "value2"), stringsAsFactors = FALSE)
table <- flextable(data = data)
table <- add_header(table, values = c("Model", "Independent variable", "p-value", "OR", "95% CI"))
table <- data.frame()
# 提取模型信息
for (i in 1:length(models)) {
  model <- models[[i]]
  model_name <- paste0("model", i)
  independent_var <- as.character(formula(model)[[2]])
  p_value <- summary(model)$coefficients[2, 4]
  odds_ratio <- exp(coef(model)[2])
  ci <- confint(model)[2, ]
  
  new_row <- data.frame(Model = model_name,
                        `Independent variable` = independent_var,
                        `p-value` = p_value,
                        OR = odds_ratio,
                        `95% CI` = paste0("[", ci[1], ", ", ci[2], "]"))
  
  table <- rbind(table, new_row)
}
#输出
write.csv(table,file = "rgsghsegs.csv")
