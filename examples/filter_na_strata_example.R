library(survival)
# create test dataset
test1 <- data.frame(case = c(0,0,1,0,1,1,0,1,NA,1),
                    x = c(0,2,1,2,1,0,1,0,0,1),
                    step_id = c(0,0,0,0,1,1,1,2,2,2))
# Remove NAs
f <- case ~ x + strata(step_id)
filter_na_strata(f, test1)

# Fit a stratified model
if(FALSE) {
  # with no NAs; necessary if using glmnet
  coxph(Surv(rep(1, length(case)), case) ~ x + strata(step_id), filter_na_strata(f, test1))
  # This differs from that
  coxph(Surv(rep(1, length(case)), case) ~ x + strata(step_id), test1)
}
