library(plyr)
library(readr)
setwd("loan data")
mydir = "csvfolder"
myfiles = list.files(pattern="*.csv", full.names=TRUE)
myfiles
# only read year 2016-2019
dat_csv = ldply(myfiles[5:19], read_csv)

all_na_row<-apply(dat_csv,1,FUN=function(x){sum(is.na(x))})

irrelevant_data = c('policy_code', 'id', 'emp_title', 'purpose', 'desc', 'title', 
                    'int_rate', 'member_id', 'acc_now_delinq', 'chargeoff_within_12_mths', 'delinq_amnt')
too_many_level<-c("earliest_cr_line","zip_code","addr_state","sub_grad","addr_state","sub_grade","emp_length")

data_leakage = c('next_pymnt_d', 'funded_amnt', 'funded_amnt_inv', 'issue_d','pymnt_plan', 
                'initial_list_status', 'out_prncp', 'out_prncp_inv', 'total_pymnt', 'total_pymnt_inv', 'total_rec_prncp',
                'total_rec_int', 'total_rec_late_fee', 'recoveries', 'collection_recovery_fee', 'last_pymnt_d', 'last_pymnt_amnt',
                'last_credit_pull_d', 'last_fico_range_high', 'last_fico_range_low', 'collections_12_mths_ex_med', 'mths_since_last_major_derog',
                'num_tl_120dpd_2m', 'num_tl_30dpd', 'hardship_flag', 'hardship_type', 'hardship_reason', 'hardship_status',
                'deferral_term', 'hardship_amount', 'hardship_start_date', 'hardship_end_date', 'payment_plan_start_date', 
                'hardship_length', 'hardship_dpd', 'hardship_loan_status', 'orig_projected_additional_accrued_interest', 
                'hardship_payoff_balance_amount', 'debt_settlement_flag', 'debt_settlement_flag_date', 'settlement_status', 
                'settlement_date', 'settlement_amount', 'settlement_percentage', 'settlement_term',
                'hardship_last_payment_amount')


lc<-dat_csv[which(all_na_row!=ncol(dat_csv)),setdiff(names(dat_csv),c(irrelevant_data,data_leakage,too_many_level))]
n<-nrow(lc)
mr<-apply(lc,2,FUN=function(x){sum(is.na(x))/n})
lc<-lc[lc$loan_status%in%c("Charged Off","Default","Fully Paid"),which(mr==0)]
lc$fico_range_avg<-(lc$fico_range_high+lc$fico_range_low)/2
lc$y<-ifelse(lc$loan_status=="Fully Paid",0,1)
lc<-lc[,setdiff(names(lc),c("loan_status","fico_range_low","fico_range_high"))]
chr<-names(lc)[sapply(lc,class)=="character"]
sapply(lc[,chr],table)
lc<-lc[lc$home_ownership!="NONE",]
lc$home_ownership[lc$home_ownership=="ANY"]<-"RENT"


write.csv(lc,"lc16_19.csv",row.names = F)





