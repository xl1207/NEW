library(broom)
library(readxl)
library(flextable)
library(officer)

# 读取Excel文件，将第一行作为列名
mydata <- read_excel("454加34填充好了.xlsx", sheet = "Sheet1", col_names = TRUE)
mydatacopy <- mydata

# 创建一个空的列表，用于存储模型
models <- list()

# 创建一个空的数据框，用于存储结果
result <- data.frame()

# 循环遍历自变量列
for (i in 2:ncol(mydatacopy)) {
  y <- unlist(mydatacopy[, 1])
  x <- unlist(mydatacopy[, i])
  
  # 进行逻辑回归
  model <- glm(y ~ x, data = mydatacopy, family = binomial(link = "logit"))
  models[[i]] <- model
  
  # 提取模型信息
  model_name <- paste0("model", i)
  independent_var <- colnames(mydatacopy)[i]
  p_value <- summary(model)$coefficients[2, 4]
  odds_ratio <- round(exp(coef(model)[2]), 4)
  ci <- round(confint(model)[2, ], 4)
  
  new_row <- data.frame(Model = model_name,
                        `Independent variable` = independent_var,
                        `p-value` = p_value,
                        OR = odds_ratio,
                        `95% CI` = paste0("[", ci[1], ", ", ci[2], "]"))
  
  result <- rbind(result, new_row)
}

# 将结果写入CSV文件
write.csv(result, file = "222.csv", row.names = FALSE)


