# exercise 1.1 
View(attenu)
attenu[is.na(attenu$station),]
attenu_cleaned <- attenu[!(is.na(attenu$station)),]
View(attenu_cleaned)
head(attenu_cleaned)
dim(attenu_cleaned)

# exercise 1.2
View(Theoph)
Theoph_2 <- Theoph
View(Theoph_2)
str(Theoph_2)
median(Theoph_2$Dose)

Theoph_2$Dose_Class <- ifelse(Theoph_2$Dose >= 4.53, "high", "low") 
head(Theoph_2)
dim(Theoph_2)

# exercise 1.3
getwd()
setwd("/Users/sabrinazhong/Desktop")
setwd("qbio_data_analysis_sabrinaz/week3_homework")
starbucks <- read.csv(file = "starbucks.csv")
View(starbucks)

is.na(starbucks)
is_row_empty <- rowSums(is.na(starbucks[,2:7])) == 0
is_row_empty
length(is_row_empty) #177
nrow(starbucks) #177
starbucks_cleaned <- starbucks[is_row_empty == TRUE,]
starbucks_cleaned

plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, xlab = "Carbohydrates (g)", ylab = "Calories", 
     main = "Carbs vs. Calories in Starbucks Drinks")

max(starbucks_cleaned$Calories) #430
starbucks_cleaned[starbucks_cleaned$Calories == 430,]
#Starbucks® Signature Hot Chocolate

plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, xlab = "Carbohydrates (g)", ylab = "Calories", 
     main = "Carbs vs. Calories in Starbucks Drinks", col=ifelse(starbucks_cleaned$Calories == 430, "red", "black"))
starbucks_cleaned$is_highest_fat <- ifelse(starbucks_cleaned$Fat == max(starbucks_cleaned$Fat), "TRUE", "FALSE")
View(starbucks_cleaned)
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, xlab = "Carbohydrates (g)", ylab = "Calories", 
     main = "Carbs vs. Calories in Starbucks Drinks", col=ifelse(starbucks_cleaned$is_highest_fat == "TRUE", "red", "black"))

#Exercise 1.4
getwd()
baseball <- read.csv(file = "Batting.csv")
View(baseball)
nrow(baseball[baseball$HR > 3,])

plot(baseball$yearID, baseball$HR, ylab = "homeruns", xlab = "year", main = "homeruns vs. year")

LA_Angels <- baseball[baseball$teamID == "LAA",]
View(LA_Angels)
plot(LA_Angels$yearID, LA_Angels$HR, ylab = "homeruns", xlab = "year", main = "homeruns vs. year for LA Angels")

ATL_PIT <- baseball[baseball$teamID == "ATL" | baseball$teamID == "PIT",]
plot(ATL_PIT$yearID, ATL_PIT$HR, ylab = "homeruns", xlab = "year", main = "homeruns vs. year for ATL and PIT",
     col=ifelse(ATL_PIT$teamID == "ATL", "red", "blue"))

#Exercise 1.5
easy_plot <- function(x, y, color_data) {
  color_med <- median(color_data)
  levels <- c(ifelse(color_data < color_med, "low", "high"))
  levels = factor(levels)
  plot(x, y, col = levels, psch = 20)
  print(color_med)
}

easy_plot(starbucks_cleaned$Calories, starbucks_cleaned$Fat, starbucks_cleaned$Sodium)

#Exercise 2.1
View(iris)
#The data set shows the sepal length, sepal width, petal length, and petal width
#of three different species of irises: setosa, versicolor, and virginica. There
#are 5 variables per data point. 

#Exercise 2.2
#The lengths and widths of the sepal and petals are continuous variables. The
#species of iris is the categorical variable. 

#Exercise 2.3
hist(iris$Sepal.Length)
hist(iris$Sepal.Width)
hist(iris$Petal.Length)
hist(iris$Sepal.Width)
#For all the histograms except the one for the petal length, the histograms have
#a bell shape to them. For the petal length, there is a very flat bell shape from
#around 3 to 7 cm. Then, below 2 cm, there is a sudden spike.  

#Exercise 2.4
mean_sepal_width <- mean(iris$Sepal.Width)
iris_copy <- iris

sepal_compare <- ifelse(iris$Sepal.Width > mean_sepal_width, "wide", "narrow")
  
iris$comparison <- sepal_compare

boxplot(iris$Sepal.Width ~ iris$comparison, xlab = "Comparison with Mean", ylab = "Sepal Width (cm)", 
        main = "Sepal Width in Comparison with the Mean Width")

#Exercise 2.5
?pairs
pairs(iris[,1:4],col=iris[,5])

#Exercise 3.1
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("TCGAbiolinks")







