#======================================StaireCase_4Point_Average=======================================
# When you feed a ".txt" file that contains the zone test data, it averages the last four point before the starecase converge
# for both negative zone and positive zone of good stereo peception.
# The average that you get "print(Average_PosZONE) & print(Average_NegZONE)" at the end is to be considered.
# Please choose this file "KZR004ZONE0.1D_FOUR_Age_ExptA_Block3_Sep_0.1_20190831T172445.txt "when prompted.
 #=====================================Code Starts Here==================================================
# Installing
#install.packages("readr")
# Loading
library(readr)

# Load data form 9th row to avoid the demografic information such as age name and participant I.D.
my_data <- read_tsv(file.choose(),col_names = TRUE,skip = 8)

#Sorting data into positive zone distance and negetive zone distances:
newdata <- my_data[order(my_data$Distance),]

#Negetive Zone
Neg_data<-subset(my_data,my_data$Distance == -1 & my_data$`Reversal?` == 1)

# Positive Zone
Pos_data<-subset(my_data,my_data$Distance == 1 & my_data$`Reversal?` == 1)

Average_PosZONE<-colMeans(Pos_data[c(9,10,11,12),11])

Average_NegZONE<-colMeans(Neg_data[c(9,10,11,12),11])

NegZONE_data<-subset(my_data,my_data$Distance == -1 )
PosZONE_data<-subset(my_data,my_data$Distance == 1 )

plot(NegZONE_data$Trial,NegZONE_data$conflict,pch = 8, col = "red",type="l", ylim = c(3,-3),xlab = "Trials",ylab = "conflict")
points(PosZONE_data$Trial,PosZONE_data$conflict, pch = 6,type="l", col = "blue")

print(Average_PosZONE)
print(Average_NegZONE)


