clinic_read_in <- read.csv("/Users/sabrinazhong/Desktop/coad_clinical_data.csv")

# Written Activity
# 1 
# A categorical variable is a category that does not relate to a value and refers to a qualitative aspect of an object. An exmaple of this is a person's sex.
# A discrete variable is a variable that is countable whole number variable. An example of this is how many days since someone's cancer diagnosis. 
# A continuous variable is a variable obtained through measuring which can be inbetween numbers. An example of this is a person's height or weight. 

# 2
colnames(clinic_read_in)
is.na(clinic_read_in$number_of_first_degree_relatives_with_cancer_diagnosis)

# Echo and I chose the "Number of 1st degree relatives with cancer diagnosis" column. 

# 3
# The variable can be collected by looking at someone's family health history or through surveys. It is a discrete variable
# because we count people as whole numbers. 

# 4
# https://www-clinicalkey-com.libproxy2.usc.edu/#!/content/playContent/1-s2.0-S095980490900344X
# This research article found that the patients who have a first degree relative diagnosed with 
# lung cancer or gastric cancer had a 4.6-fold higher chance of having an early onset of non-small
# cell lung cancer (NSCLC). They believe the reason for this is due to genetic factors or having
# similar environmental exposures with their relatives. 

# https://ashpublications.org/blood/article/134/12/960/374916/Analysis-of-153-115-patients-with-hematological
# Those with a first degree relative that has blood cancer also had a higher risk of having blood
# cancer. This could be due to how blood cancer is caused by heritable DNA changes. 

# 5
is.na(clinic_read_in$race_list) 
# We chose the ethnicity variable. This data can be collected through a simple questionnaire. This
# is a categorical variable because it is a qualitative characteristic of a patient.

# 6
# Minorities might have more first degree relatives with cancer because they might be more likely to
# be in environments that can induce stress to their health or pose hazards to their health.
# Patients with more first degree relatives with CRC will have a lower survival rate.
# Minority patients will have a lower survival rate. 

# 7
# Not sure if I created the graph correctly, but it seems as though the graph I created illustrates
# the opposite of the research article findings. The patients with no first degree relatives with 
# cancer have a lower survival probability than those with the first degree relatives with cancer. 

# In the graph showing the race list. There does seem to be a difference in survival for each
# race. The minorities seem to have a lower probability of survival except for the Asians. 

par(mar=c(5,4,3,1))
boxplot(clinic_read_in$number_of_first_degree_relatives_with_cancer_diagnosis ~ clinic_read_in$race_list, 
        xlab = "Race", ylab = "Number of 1st Degree Relatives with Cancer", main = "Number of 1st Degree 
        Relatives with Cancer vs. Race", las = 2)


install.packages("survival")
library(survival)
install.packages("survminer")
library(survminer)
clinic_read_in$death_event <- ifelse(clinic_read_in$vital_status == "Alive", 0, 1)
# initializing the survival object
surv_object <- Surv(time = clinic_read_in$days_to_death, 
                    event = clinic_read_in$death_event)

# KM plot for # of 1st degree relatives
relatives_fit <- surv_fit( surv_object ~ clinic_read_in$number_of_first_degree_relatives_with_cancer_diagnosis, data = clinic_read_in )

#formatting of the plot
survplot = ggsurvplot(relatives_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "bottom")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=10), # increase font sizes
        axis.text = element_text(size=9),
        legend.title = element_text(size=10),
        legend.text = element_text(size=9))
p

# KM plot for Race (sorry didn't realize it would be the same as last hw)
# creating object for race fit
race_fit <- surv_fit( surv_object ~ clinic_read_in$race_list, data = clinic )

# formatting
survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "bottom")
# formatting
p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=10), # increase font sizes
        axis.text = element_text(size=9),
        legend.title = element_text(size=10),
        legend.text = element_text(size=9))
p













