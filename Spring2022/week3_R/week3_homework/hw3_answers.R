# exercise 1.1
attenu
noStation <- !is.na(attenu$station)
attenu_cleaned <-attenu[noStation,]
print(head(attenu_cleaned))
print(dim(attenu_cleaned))

#exercise 1.2
Theoph_2 <- Theoph
str(Theoph_2)
median(Theoph_2[,3])
Theoph_2$Dose_Class <- ifelse(Theoph_2[,3] >= 4.54, "high", "low")
print(head(Theoph_2))
print(dim(Theoph_2))

#exercise 1.3 
setwd("~/Programming_Projects/QBIOPublicDataAnalysis/qbiodataanalysisleon/Spring2022/week3_R/week3_homework")
starbucks <- read.csv("starbucks.csv")
na_dataframe <- is.na(starbucks)
is_row_empty <- ifelse(rowSums(na_dataframe[,]) == 0,FALSE, TRUE )
print(nrow(starbucks) == length(is_row_empty))
starbucks_cleaned <- starbucks[!is_row_empty,]
plot(x = starbucks_cleaned$Calories, 
     y = starbucks_cleaned$Carb,
     main = "Starbucks Drinks Calories vs Carbs",
     xlab = "Calories",
     ylab = "Carbs (grams)")
#calories and carbs seem correlated
highest_calorie <- max(starbucks_cleaned$Calories)
rowIndex_highest_calorie <- starbucks_cleaned[starbucks_cleaned$Calories == highest_calorie,]
dim(rowIndex_highest_calorie)
rowIndex_highest_calorie$Drink
starbucks_cleaned$is_highest_fat <- (starbucks_cleaned$Calories == highest_calorie)
plot(x = starbucks_cleaned$Calories, 
     y = starbucks_cleaned$Carb,
     col = factor(starbucks_cleaned$is_highest_fat),
     main = "Starbucks Drinks Calories vs Carbs",
     xlab = "Calories",
     ylab = "Carbs (grams)")


#exercise 1.4
baseball <- read.csv("Batting.csv")
head(baseball)
scoredAtLeastThree <- length(baseball[baseball$HR > 3,])
scoredAtLeastThree
# 
# 
# plot_baseball = function(teamName, dataframe) {
#   agg <- aggregate(x = dataframe[, 12],
#             by = list(dataframe$yearID),
#             FUN = count)
#   plot(x = agg$Group.1, y = agg$HR)
# }
# 
# plot_baseball("All Teams", baseball)
# 
# LA_Angels <- baseball[baseball$teamID == "LAA", ]
# plot_baseball("LA Ange;s", LA_Angels)
# 
# two_teams <- baseball[baseball$teamID == "PIT" | baseball$teamID == "ATL",]
# 

plot_baseball = function(teamName, dataframe) {
  plot(x = dataframe$yearID, y = dataframe$HR, col = factor(dataframe$teamID),
       xlab = "Year", ylab = "Num Home Runs By A Player", main = paste("Home Runs Over Time Across", teamName))
}

plot_baseball("All Teams", baseball)
LA_Angels <- baseball[baseball$teamID == "LAA", ]
plot_baseball("LA Angels", LA_Angels)

Pit_or_Atl <- baseball[baseball$teamID == "PIT" | baseball$teamID == "ATL", ]
plot_baseball("The PIT And ATL Teams", Pit_or_Atl)


#exercise 1.5
easy_plot = function(x, y, color_data) {
  median_color <- median(color_data)
  levels <- ifelse(color_data < median_color, "low","high")
  levels = factor(levels)
  print(median_color)
  
  plot(x = x, y = y, col = levels, pch=20)
  print(cor.test(x, y))
}

easy_plot(c(1, 2, 3), c(4, 5, 6), c(1, 2, 3))
#each datapoint is a different color
easy_plot(starbucks_cleaned$Calories, starbucks_cleaned$Carb, starbucks_cleaned$Sodium)
#yes, there are significant correlations


#Exercise 2.1
iris
#This data gives measurements for flower characteristics of various species of irises
dim(iris)
#150 observations
#5 features per observation

#Exercise 2.2
library(tidyverse)
glimpse(iris)
#sepal length, sepal width, petal length, and petal width are continous. Species is categorical. The types are dbl, dbl, dbl, dbl, fct, respectively
hist(iris$Sepal.Length)
hist(iris$Sepal.Width)
hist(iris$Petal.Length)
hist(iris$Petal.Width)
#the sepal distributions are appropximately normal, while the petal distribution includes a large amount of observations in the lower end
mean_sepal_width <- mean(iris$Sepal.Width)
iris_copy <- iris
narrow_sepal <- ifelse(iris_copy$Sepal.Width < mean_sepal_width, TRUE, FALSE)
iris_copy$narrow_sepal <- narrow_sepal
boxplot(iris_copy$Sepal.Width ~ iris_copy$narrow_sepal)

#Exercise 2.5
#based on the plots, setosa seems to be the most unique, while versicolor and virginica are more similar
iris
?pairs
args(pairs)
pairs(iris_copy[, -5:-6])
