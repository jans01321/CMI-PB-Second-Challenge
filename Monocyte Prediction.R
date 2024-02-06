setwd('/Multi-Omics data challenge/Code')
###### monocyte1 #######
# read data and identify training and testing data
data_imput = read.csv('/Multi-Omics data challenge/Code/imputed_data_72.csv')
data_imput = data_imput[,-1]
train_id = which(data_imput$train_set==1)
predictor_mon1 = read.table('/Multi-Omics data challenge/Code/predictors.txt')
predictor_mon1 = as.vector(predictor_mon1[,1])
predictor_mon1 = predictor_mon1[-1]
X = data_imput[train_id,predictor_mon1]
Y = data_imput[train_id,'Monocytes_1']
test_X = data_imput[-train_id,predictor_mon1]
test_Y = data_imput[-train_id,'Monocytes_1']

# fit elastic net model with lambda tuned by cross-validation
library(glmnet)
fit.en = cv.glmnet(data.matrix(X), as.matrix(Y), family="gaussian", alpha=0.5)
bestlambda<-fit.en$lambda.min
test_result <- predict(fit.en, s=bestlambda, newx = data.matrix(test_X)) # coef were estimated from the entire training set
true_rank = rank(test_Y)
test_rank = rank(test_result)

cor(true_rank, test_rank, method = 'spearman')

# obtain predicted rank from test data set
# read data
cell = read.csv('pbmc_cell_frequency_test.csv', row.names = 1)
cell=as.data.frame(t(cell))
cell$ID = row.names(cell)
gene = read.csv('pbmc_gene_expression_test.csv', row.names = 1)
gene = as.data.frame(t(gene))
gene$ID = row.names(gene)
anti = read.csv('plasma_antibody_levels_test.csv', row.names = 1)
anti = as.data.frame(t(anti))
anti$ID = row.names(anti)
cyto = read.csv('plasma_cytokine_concentrations_test.csv', row.names = 1)
cyto = as.data.frame(t(cyto))
cyto$ID = row.names(cyto)

merged_df1 <-  merge(cell, gene, by = "ID", all=T)
merged_df2 <-  merge(merged_df1, anti, by = "ID", all=T)
final_merged_df <- merge(merged_df2, cyto, by = "ID", all=T)

final_merged_df$ID = gsub("X", "", final_merged_df$ID)
saveRDS(final_merged_df, 'merged_df.rds')

# load imputed data and rename
final_merged_df = read.csv('imputed_test_63.csv')

pred_X = final_merged_df[,predictor_mon1]
pred_result <- predict(fit.en, s=bestlambda, newx = data.matrix(pred_X)) # coef were estimated from the entire training set
pred_rank_mon1 = rank(pred_result)

###### fold change ######
library(dplyr)
data_imput <- data_imput %>% mutate(change = Monocytes_1/Monocytes_0)
predictor = read.csv('predictors_change.txt')
predictor_change = as.vector(unlist(predictor))
X = data_imput[train_id, predictor_change]
Y = data_imput[train_id, 'change']
test_X = data_imput[-train_id, predictor_change]
test_Y = data_imput[-train_id,'change']
fit.en = cv.glmnet(data.matrix(X), as.matrix(Y), family="gaussian", alpha=0.5)
bestlambda<-fit.en$lambda.min
test_result <- predict(fit.en, s=bestlambda, newx = data.matrix(test_X)) # coef were estimated from the entire training set
true_rank = rank(test_Y)
test_rank = rank(test_result)

cor(true_rank, test_rank, method = 'spearman')

pred_X = final_merged_df[,predictor_change]
pred_result <- predict(fit.en, s=bestlambda, newx = data.matrix(pred_X)) # coef were estimated from the entire training set
pred_rank_change = rank(pred_result)

pred_result = data.frame('ID' = final_merged_df$proc_test...1., 'Rank for Monocytes1' = pred_rank_mon1,
                         'Rank for fold change' = pred_rank_change)
write.csv(pred_result, 'pred_result.csv', row.names = F)
