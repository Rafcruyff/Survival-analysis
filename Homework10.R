library(survival)
?cgd
data(cgd)

cox_model <- coxph(Surv(time,death), data = cgd)
