##############################################
# Oncology Risk Modeling and Export Pipeline
# Data: TCGA MAF + Clinical (Colorectal, 2025)
# All outputs will be saved in:
#   G:/Oncology_Algorithm/cohortMAF.2025-05-25.maf/
##############################################

# 0. Setup
library(data.table)
library(dplyr)
library(ggplot2)
library(pROC)
library(survival)
library(survminer)
library(randomForest)

outdir <- "G:/Oncology_Algorithm/cohortMAF.2025-05-25.maf/"
dir.create(outdir, showWarnings = FALSE)

# 1. LOAD DATA FILES --------------------------
maf_file <- "G:/Oncology_Algorithm/cohortMAF.2025-05-25.maf/cohortMAF.2025-05-25.maf"
clinical_file <- "G:/Oncology_Algorithm/clinical.cohort.2025-05-25/clinical.tsv"

maf <- fread(maf_file)
clinical <- fread(clinical_file)

# 2. FEATURE AGGREGATION ----------------------
drivers <- c("KRAS", "TP53", "APC", "SMAD4")

flags <- dcast(
  maf[Hugo_Symbol %in% drivers &
        Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins")],
  Tumor_Sample_Barcode ~ Hugo_Symbol,
  fun.aggregate = function(x) as.integer(length(x) > 0),
  value.var = "Hugo_Symbol"
)

tmb <- maf[
  Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins"),
  .N,
  by = Tumor_Sample_Barcode
]
setnames(tmb, "N", "Mutation_Count")

flags[, Sample_ID := substr(Tumor_Sample_Barcode, 1, 12)]
tmb[, Sample_ID := substr(Tumor_Sample_Barcode, 1, 12)]

feature_table <- merge(flags, tmb, by = c("Tumor_Sample_Barcode", "Sample_ID"), all = TRUE)
feature_table <- merge(feature_table, clinical, by.x = "Sample_ID", by.y = "cases.submitter_id", all.x = TRUE)

# 3. SELECT, CLEAN, FORMAT --------------------
feature_table$Age <- suppressWarnings(as.numeric(as.character(feature_table$demographic.age_at_index)))

model_data <- feature_table[, .(
  KRAS, TP53, APC, SMAD4, Mutation_Count,
  Age,
  Stage = diagnoses.ajcc_pathologic_stage,
  Outcome = demographic.vital_status,
  FollowUp = diagnoses.days_to_last_follow_up
)]

model_data$Stage <- as.factor(model_data$Stage)
model_data$Outcome <- as.factor(model_data$Outcome)
model_data <- na.omit(model_data)

# 4. TRAIN/TEST SPLIT -------------------------
set.seed(1234)
train_idx <- sample(seq_len(nrow(model_data)), size = 0.7 * nrow(model_data))
train <- model_data[train_idx, ]
test <- model_data[-train_idx, ]

# 5. LOGISTIC REGRESSION ----------------------
logit_model <- glm(Outcome ~ KRAS + TP53 + APC + SMAD4 + Mutation_Count + Age + Stage,
                   data = train, family = binomial)
summary(logit_model)

# 6. ROC, CONFUSION, EXPORT -------------------
pred_probs <- predict(logit_model, newdata = test, type = "response")
roc_obj <- roc(test$Outcome, pred_probs)
auc_val <- auc(roc_obj)
png(file.path(outdir, "ROC_logistic.png"), width=600, height=400)
plot(roc_obj, main = "ROC Curve (Logistic Regression)")
dev.off()

pred_class <- ifelse(pred_probs > 0.5, "Dead", "Alive")
conf_mat <- table(Predicted = pred_class, Actual = test$Outcome)
write.csv(conf_mat, file.path(outdir, "ConfusionMatrix_logistic.csv"))

# 7. RANDOM FOREST ----------------------------
rf_model <- randomForest(Outcome ~ KRAS + TP53 + APC + SMAD4 + Mutation_Count + Age + Stage,
                         data = train, ntree = 500, importance = TRUE)
rf_probs <- predict(rf_model, newdata = test, type = "prob")[,2]
rf_roc <- roc(test$Outcome, rf_probs)
auc_rf <- auc(rf_roc)
png(file.path(outdir, "ROC_randomforest.png"), width=600, height=400)
plot(rf_roc, main = "ROC Curve (Random Forest)")
dev.off()
png(file.path(outdir, "RF_VariableImportance.png"), width=600, height=400)
varImpPlot(rf_model, main = "Random Forest Variable Importance")
dev.off()

# 8. RISK STRATIFICATION ----------------------
test$risk_score <- pred_probs
test$risk_group <- cut(test$risk_score, breaks = quantile(test$risk_score, probs = c(0,0.33,0.66,1)),
                       labels = c("Low", "Medium", "High"), include.lowest = TRUE)

# 9. SURVIVAL ANALYSIS ------------------------
test_surv <- test[
  !is.na(FollowUp) & FollowUp != "" &
    !is.na(Outcome) & Outcome != "" &
    !is.na(risk_group) &
    !is.na(as.numeric(FollowUp))
]
test_surv$FollowUp <- as.numeric(test_surv$FollowUp)
fit <- survfit(Surv(time = FollowUp, event = ifelse(Outcome == "Dead", 1, 0)) ~ risk_group, data = test_surv)
kmplot <- ggsurvplot(fit, data = test_surv, risk.table = TRUE, title = "Kaplan-Meier by Risk Group")
ggsave(filename = file.path(outdir, "KM_risk_stratification.png"), plot = kmplot$plot)
ggsave(filename = file.path(outdir, "KM_risk_stratification_risktable.png"), plot = kmplot$table)

cat("Events (deaths):", sum(ifelse(test_surv$Outcome == "Dead", 1, 0)), "\n")
cat("Median follow-up (days):", median(test_surv$FollowUp, na.rm = TRUE), "\n")

# 10. EXPORT MAJOR TABLES ---------------------
exp_coef <- exp(cbind(OR = coef(logit_model), confint(logit_model)))
write.csv(exp_coef, file.path(outdir, "Logistic_ORs.csv"))
write.csv(data.frame(variable=names(coef(logit_model)), coef=coef(logit_model)),
          file.path(outdir, "Model_coefficients_for_EMR.csv"))

table1 <- train %>%
  group_by(Outcome) %>%
  summarise(
    n = n(),
    Age_mean = mean(Age, na.rm=TRUE),
    Mutation_Count_mean = mean(Mutation_Count, na.rm=TRUE),
    KRAS_mut = mean(KRAS, na.rm=TRUE),
    TP53_mut = mean(TP53, na.rm=TRUE),
    APC_mut = mean(APC, na.rm=TRUE),
    SMAD4_mut = mean(SMAD4, na.rm=TRUE)
  )
write.csv(table1, file.path(outdir, "Table1_Baseline.csv"))
fwrite(model_data, file.path(outdir, "model_data_cleaned.csv"))

##############################################
# END OF PIPELINE
##############################################

