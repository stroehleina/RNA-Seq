# Created by Andreas J Stroehlein, The University of Melbourne, 2015
# No package dependencies

##FUNCTIONS####
# function that rounds up to next highest ten
roundUp <- function(x,to=10)
{
to*(x%/%to + as.logical(x%%to))
}
##############


#START
#read in table of TPM values (attached) from RSEM output
trans_table <- read.table("clean_transcription_table.final", row.names=1, header=FALSE)

#read in kinome data table (attached). Only needed when interested in pathway information for your genes
#master <- read.table("master_Suppl_432_29_09.cleanIDs", row.names=1, header=FALSE, na.strings="", sep="\t")
#master[is.na(master[,15]),15] <- 0

#name columns, not needed but helps to check if you are dealing with lots of columns
colnames(trans_table) <- c("gene"    ,"egg_cnt" ,  "egg_tpm" ,  "L1_cnt" ,   "L1_tpm"  ,  "L2a_cnt" ,  "L2a_tpm" ,  "L2b_cnt" ,  "L2b_tpm"  , "L3_cnt" ,   "L3_tpm"  ,  "L4Fe_cnt" , "L4Fe_tpm" , "L4Ma_cnt" , "L4Ma_tpm" ,"AduFe_cnt", "AduFe_tpm" ,"AduMa_cnt" ,"AduMa_tpm" ,"Group")

#select columns with tpm values only
tpm_table <- trans_table[, c(seq(3,7,by=2),seq(11,19,by=2))]

#create distance matrix Euclidean
dist_tpm <- dist(as.matrix(tpm_table[,c(1:8)]), method="euclidean")

#run Ward clustering
clust_tpm <- hclust(dist_tpm, method="ward.D2")

#select k =15. This cuts the clustering dendogram at a certain height where number of groups = k = 15.
#Rather arbitrary value. Depends on how fine you want the groups to resolve
cluster <- cutree(clust_tpm, k=15)

#add cluster numbers as extra column to table
tpm_table <- cbind(tpm_table, cluster)

# I have used this section to add additional KEGG pathway information in the last columns of the table.
# Might be useful if you want to look at pathway enrichment in certain transcription profiles
#tpm_anno <-  trans_table[, c(seq(3,7,by=2),seq(11,19,by=2),20)]
#tpm_anno <- cbind(tpm_anno, cluster)
#pathways <- master[c(25,26,27)]
#tpm_anno <- cbind(tpm_anno, pathways)
#colnames(tpm_anno)[c(11:13)] <- c("pw2","pw3","pw4")

#create pdf
pdf(file="check_17_05.pdf")

#print clusters in 4x4 format, tweak margins
par(mfrow=c(4,4), mai=c(0.4,0.2, 0.1, 0.1), oma=c(0.5,1,0.5,0.5))

#counter variable
c <- 0
# roll your own, define the order of clusters in which you want the graphs to be output in the pdf
for(k in unique(tpm_table$cluster)[c(4,8,7,3,5,6,12,13,1,2,9,10,11,14,15)]){
c <- c+1
#calculate standard deviation
sd_trans <- apply(tpm_table[tpm_table$cluster == k,],2,sd)

#if cluster is in the first column set label to "TPM"
if(k %% 4 ==1){
lab <- "TPM"
}else{
lab <- ""
}

#PLOT
#plot graph. set y axis range to max value rounded up to next 10 (e.g., 432 becomes 440 as upper limit). don't print labels or axes
plot(0,0,xlim = c(0,10),ylim = c(0,roundUp(max(range(tpm_table[tpm_table$cluster == k,])))),type = "n", ylab = "", yaxt="n", xaxt="n",xlab="")
#add axes
axis(2, cex.axis=0.8)
axis(1,labels=c("Egg","L1","L2","L3","L4","Adult"), at=seq(0,10, by=2), las=2, cex.axis=0.8)
#add x axis description
title(xlab=paste("Cluster ",c," (n = ",nrow(tpm_table[tpm_table$cluster == k,]),")", sep=""), cex.lab = 0.8, line = 2.1)
#FIXME this doesn't seem to work properly, "TPM" is either printed in the middle of the pdf (using mtext) or not at all (using title)
title(ylab=lab)
#mtext(lab, side=2, cex = 0.5, outer=TRUE)

#select all rows with current cluster number and get values for egg, L1, L2, L3. For each gene plot line from position 0 to 6 in light gray
for(i in row.names(tpm_table[tpm_table$cluster == k,])){lines(seq(0,6, by=2),tpm_table[i,c(1,2,3,4)], col="gray90")}
#select all rows with current cluster number and get values for L3, L4 female, adult female. plot from position 6 to 10 in light red
for(i in row.names(tpm_table[tpm_table$cluster == k,])){lines(seq(6,10, by=2),tpm_table[i,c(4,5,7)], col="lightpink")}
#select all rows with current cluster number and get values for L3, L4 male, adult male. plot from position 6 to 10 in light blue
for(i in row.names(tpm_table[tpm_table$cluster == k,])){lines(seq(6,10, by=2),tpm_table[i,c(4,6,8)], col="lightblue")}

#TRENDLINES
#Calculate and print trendlines using Lowess regression method for female, male and EGG,L1-L3

#get only female columns
female_ts <- rbind(cbind(tpm_table[tpm_table$cluster == k,4],"6"),cbind(tpm_table[tpm_table$cluster == k,5],"8"),cbind(tpm_table[tpm_table$cluster == k,7], "10"))
#Lowess function
female_lo <- lowess(female_ts[,2], female_ts[,1])
#print trendline
lines(unique(female_lo$x), unique(female_lo$y), col="red", lwd=1)
#print upper and lower dashed lines representing +-standard deviation female 
upperfe <- unique(female_lo$y) + sd_trans[c(4,5,7)]
lines(seq(6,10, by=2), upperfe, col="red" , lwd=1, lty=2)
lowerfe <- unique(female_lo$y) - sd_trans[c(4,5,7)]
lines(seq(6,10, by=2), lowerfe, col="red" , lwd=1, lty=2)

#get only male columns
male_ts <- rbind(cbind(tpm_table[tpm_table$cluster == k,4],"6"),cbind(tpm_table[tpm_table$cluster == k,6],"8"),cbind(tpm_table[tpm_table$cluster == k,8], "10"))
#Lowess function
male_lo <- lowess(male_ts[,2], male_ts[,1])
#print trendline
lines(unique(male_lo$x), unique(male_lo$y), col="blue", lwd=1)
#print upper and lower dashed lines representing +-standard deviation male
upperma <- unique(male_lo$y) + sd_trans[c(4,6,8)]
lines(seq(6,10, by=2), upperma, col="blue" , lwd=1, lty=2)
lowerma <- unique(male_lo$y) - sd_trans[c(4,6,8)]
lines(seq(6,10, by=2), lowerma, col="blue" , lwd=1, lty=2)

#get only Egg to l3 columns
l1l3_ts <- rbind(cbind(tpm_table[tpm_table$cluster == k,1],"0"),cbind(tpm_table[tpm_table$cluster == k,2],"2"),cbind(tpm_table[tpm_table$cluster == k,3], "4"),cbind(tpm_table[tpm_table$cluster == k,4],"6"))
#Lowess function
l1l3_lo <- lowess(l1l3_ts[,2], l1l3_ts[,1])
#print trendline
lines(unique(l1l3_lo$x), unique(l1l3_lo$y), col="black", lwd=1)
#print upper and lower dashed lines representing +-standard deviation egg-l3
upperl1l3 <- unique(l1l3_lo$y) + sd_trans[c(1,2,3,4)]
lines(seq(0,6, by=2), upperl1l3, col="black" , lwd=1, lty=2)
lowerl1l3 <- unique(l1l3_lo$y) - sd_trans[c(1,2,3,4)]
lines(seq(0,6, by=2), lowerl1l3, col="black" , lwd=1, lty=2)

}

dev.off()
#END