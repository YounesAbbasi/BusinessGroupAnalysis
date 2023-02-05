#"This R Script for "The structure and formation of business group:Evidence from Korean chaebols" article for Iran

setwd ("D:/R_Project/Business_Group_project")

#######1)Ultimate Cashflow Right, postition and loops
###### Importing ownership matrix and f vector
####importing malekiat matrix that is output of python network code

malekiat_Data_Frame = read.csv("D:/R_Project/Business_Group_project/processed/malekiat_matrix_of_family.csv", sep = ",")

####Deleting first column
malekiat_Data_Frame[1]=NULL

####converting data to numeric for calculations
malekiat_Data_Frame[,1:ncol(malekiat_Data_Frame)] = lapply(malekiat_Data_Frame[,1:ncol(malekiat_Data_Frame)], as.numeric)

####converting dataframe to matrix for calculations
malekiat_matrix = as.matrix(malekiat_Data_Frame)

#### importing head's direct share of any firms
direct_share_of_family = read.csv("f_Vector_for_family.csv")

#### delete columns that aren't in use
direct_share_of_family1=subset(direct_share_of_family,select = -c(1,2,4))

#### converting direct share to matrix
direct_share_of_family_vector = as.matrix(direct_share_of_family1)

###### Generating ultimate ownership
Ultimate_ownership=t(t(direct_share_of_family_vector) %*% solve(diag(ncol(malekiat_matrix))-malekiat_matrix))

#### list of firms that head has ultimate ownership in it
lsit_of_ownership1 = which(Ultimate_ownership %in% Ultimate_ownership[Ultimate_ownership>0])

#### taking the national code of tims in list of ownership

list_of_ownership_names = c()
for (i in lsit_of_ownership1) (list_of_ownership_names = c(list_of_ownership_names,direct_share_of_family[i,1]))

#### taking output from U and national code

list_of_ownership_and_u_value = cbind(list_of_ownership_names,Ultimate_ownership[Ultimate_ownership>0])

write.csv(list_of_ownership_and_u_value, file="D:/R_Project/Business_Group_project/U_Vector_and_national_Code.csv")

###### position of firm

#### positen function

position_i = function (i){
  u_i = Ultimate_ownership[i,1]
  d_i = matrix(0,nrow(direct_share_of_family_vector), 1)
  d_i[i] = 1
  Position_i_is = (1/u_i)*t(direct_share_of_family_vector) %*% (solve(diag(ncol(malekiat_matrix))-malekiat_matrix) %*% solve(diag(ncol(malekiat_matrix))-malekiat_matrix))%*%d_i
  return(Position_i_is)
}

star_time = Sys.time()

postition_of_ultimated = lapply(1:nrow(direct_share_of_family_vector), position_i)

end_time = Sys.time()


end_time - star_time
###### Firms in loop

#### loop function

library (matrixcalc)
loop_i = function(i){
  d_i = matrix(0,nrow(direct_share_of_family_vector),1)
  d_i[i] = 1
  all_loops=c()
  for (n in 1:nrow(direct_share_of_family_vector)) {
    if (t(d_i) %*% matrix.power(malekiat_matrix,n)%*%d_i>0) {
      all_loops=c(all_loops,n)
      break
    }
  }
  return(all_loops)}


loop_i_for_any_firm = lapply(1:nrow(direct_share_of_family_vector), loop_i)

#######2)Control Right and Centrality

##############function for set of firms under threshold

function_of_C_T = function(Threshold) {
  
  s_1=c()
  for (i in 1:nrow(direct_share_of_family_vector)) {
    if (direct_share_of_family_vector[i]>Threshold) {
      s_1 = unique(c(s_1,i))
    }
  }

##################Repeat to take set of firms that head has voting right in them subject to threshold
  s_2 = c()
  s_j_i = c()
  repeat{
    for (i in 1:nrow(direct_share_of_family_vector)) {if (i < nrow(direct_share_of_family_vector)){
      for (j in s_1){
        if (direct_share_of_family_vector[i]+sum(c(s_j_i,malekiat_matrix[c(j),i]))>Threshold) {
          s_2 = unique(c(s_2,i))}
      }
    }}
  
    if (all(length(s_2) == length(s_1) && all(s_2 == s_1)) == TRUE){
      break}
    s_1=s_2
    }
    return(s_1)
}




###### The critical control threshold is the highest control threshold that is consistent with family control of firm i. (for past)
critical_control_threshold_function = function (steps_of_sequence){
sequence_of_Threshold = as.array(seq(0,1,by=steps_of_sequence))
out = vector("list", length(sequence_of_Threshold))
cci = vector("list", length(sequence_of_Threshold))
  
  
start_time = Sys.time()
for (i in seq_along(sequence_of_Threshold)) {
  if (length(function_of_C_T(sequence_of_Threshold[i]))<=length(function_of_C_T(0)) && length(function_of_C_T(sequence_of_Threshold[i]))>1){
    out[[i]]=function_of_C_T(sequence_of_Threshold[i])}}
end_time = Sys.time()
  
end_time-start_time

for (i in seq_along(sequence_of_Threshold)) {skip_to_next3=FALSE
tryCatch(if(length(setdiff(as.matrix(data.frame(out[i])[,1]),as.matrix(data.frame(out[i+1])[,1])))>0)
{cci[[i]]=setdiff(as.matrix(data.frame(out[i])[,1]),as.matrix(data.frame(out[i+1])[,1]))},error=function(e){if(skip_to_next3){next}})}

max_Treshol = data.frame(colnames(c("max_T_for_any_firm","index_of_firm")))
for (i in seq_along(sequence_of_Threshold)){
  if (typeof(as.array(cci[i])[[1]])=="NULL") {
    next
  }
  
  w = data.frame( "max_T_for_any_firm"=sequence_of_Threshold[i+1],"index_of_firm"=as.array(cci[i])[[1]])
  max_Treshol = rbind.data.frame(max_Treshol,w)}
max_Treshol

i = length(seq_along(sequence_of_Threshold))
while(i>0){
  i = i-1
  
  if(typeof(as.array(out[i])[[1]]) !="NULL") {
    last_items_in_out = data.frame(list(a = c(out[i])))
    break}}

names(last_items_in_out)=c("index_of_firm")
last_items_in_out[,"max_T_for_any_firm"]=1
max_Treshol=rbind.data.frame(max_Treshol,last_items_in_out)
rownames(max_Treshol)=1:nrow(max_Treshol)
return(max_Treshol)
}

#########################################################################################
################## Voting Right function for firms in C_T

VR_i_T_function = function(i,Threshold) {
  s_j_i_1 = c()
  for (j in function_of_C_T(Threshold)){
    s_j_i_1 = c(s_j_i_1,malekiat_matrix[c(j),i])
    VR_i_T = direct_share_of_family_vector[i]+sum(s_j_i_1)}
  return (VR_i_T)
}

voting_right_function = function(Threshold) {
  VR_i_for_any_firm_in_C_T = mapply(VR_i_T_function, function_of_C_T(Threshold) ,c(Threshold))
  return (VR_i_for_any_firm_in_C_T)
}

index_and_VR_function = function(Threshold) {
  Index_and_VR = cbind(function_of_C_T(Threshold), voting_right_function(Threshold))
  print ("first column is index of f_vector and second is voting right")
  return (Index_and_VR)
}

###########################################################################

######### Centrality of firm i
#To calculate centrality of firm i functions was repeated.

#####################################################

critical_control_threshold_function_minus_firm_i = function(frim_i,sequence_by){
  
  malekiat_matrix_zero_i = malekiat_matrix
  
  malekiat_matrix_zero_i[frim_i,] = 0
  
  function_of_C_T_del_i = function(Threshold) {
    
    s_1_del_i = c()
    for (i in 1:nrow(direct_share_of_family_vector)){
      if (direct_share_of_family_vector[i]>Threshold){
        s_1_del_i = unique(c(s_1_del_i,i))
      }
    }
    
    ################# Repeat to take set of firms that head has voting right in them subject to threshold
    s_2_del_i = c()
    s_j_i = c()
    
    repeat{
      for (i in 1:nrow(direct_share_of_family_vector)) {if (i < nrow(direct_share_of_family_vector) ) {
        for (j in s_1_del_i){
          if (direct_share_of_family_vector[i]+sum(c(s_j_i,malekiat_matrix_zero_i[c(j),i]))>Threshold){
            s_2_del_i = unique(c(s_2_del_i,i))}
        }
      }}
      
      if(all(length(s_2_del_i) == length(s_1_del_i) && all(s_2_del_i == s_1_del_i)) == TRUE) {
        break}
      s_1_del_i = s_2_del_i
    }
    return(s_1_del_i)
  }
  ###### The critical control threshold is the highest control threshold that is consistent with family control of firm i. (for past)
sequence_of_Threshold = as.array(seq(0,1,by=sequence_by))
out_del_i = vector("list", length(sequence_of_Threshold))
cci_del_i = vector("list", length(sequence_of_Threshold))
  
start_time = Sys.time()
  
for (i in seq_along(sequence_of_Threshold)){
  if(length(function_of_C_T_del_i(sequence_of_Threshold[i]))<=length(function_of_C_T_del_i(0))&&length(function_of_C_T_del_i(sequence_of_Threshold[i]))>1){
    out_del_i[[i]]=function_of_C_T_del_i(sequence_of_Threshold[i])}}
  
end_time = Sys.time()
  
end_time-start_time
  
for (i in seq_along(sequence_of_Threshold)) {skip_to_next2=FALSE
  tryCatch(if(length(setdiff(as.matrix(data.frame(out_del_i[i])[,1]),as.matrix(data.frame(out_del_i[i+1])[,1])))>0)
{cci_del_i[[i]]=setdiff(as.matrix (data.frame(out_del_i[i])[,1]),as.matrix(data.frame(out_del_i[i+1])[,1]))}, error=function(e){if(skip_to_next2){nex}})}

max_Treshol_del_i = data.frame(colnames(c("max_T_for_any_firm","index_of_firm")))
for (i in seq_along(sequence_of_Threshold)){
  if(typeof(as.array(cci_del_i[i])[[1]]) =="NULL"){
    next
  }
  w = data.frame( "max_T_for_any_firm"=sequence_of_Threshold[i+1], "index_of_firm"=as.array(cci_del_i[i])[[1]])
  max_Treshol_del_i = rbind.data.frame(max_Treshol_del_i,w)}


i = length(seq_along(sequence_of_Threshold))
while (i>0) {
  i = i-1
  
  if(typeof(as.array(out_del_i[i])[[1]]) !="NULL"){
    last_items_in_out = data.frame(list(a = c(out_del_i[i])))
    break}}

names(last_items_in_out)=c("index_of_firm")
last_items_in_out[,"max_T_for_any_firm"]=1
max_Treshol_del_i=rbind.data.frame(max_Treshol_del_i,last_items_in_out)
rownames(max_Treshol_del_i)=1:nrow(max_Treshol_del_i)
#max_Treshol_del_i
return(max_Treshol_del_i)
}

central_function_for_firm_i = function(firm_i,steps_of_sequence){
  central_i = ((sum(critical_control_threshold_function(steps_of_sequence)[,1])-critical_control_threshold_function(steps_of_sequence)[match(firm_i,critical_control_threshold_function(steps_of_sequence)[,2]),1])-(sum(critical_control_threshold_function_minus_firm_i(firm_i,steps_of_sequence)[,1])-critical_control_threshold_function_minus_firm_i(firm_i,steps_of_sequence)[match(firm_i,critical_control_threshold_function_minus_firm_i(firm_i,steps_of_sequence)[,2]),1]))/(nrow(direct_share_of_family_vector)-1)
  return(central_i)}

centrality_for_all_firms_with_index_function = function (steps_of_sequence) {
  centerality_of_all_firms = mapply(central_function_for_firm_i,1:nrow(direct_share_of_family_vector), c(steps_of_sequence))
  centerality_of_all_firms = cbind(1:nrow(direct_share_of_family_vector), centerality_of_all_firms)
  return (centerality_of_all_firms)
}

###########Functions in use

#1 
Ultimate_ownership
#2
postition_of_ultimated
#3
loop_i_for_any_firm
#4
function_of_C_T(0.2) #Variable of this function is Threshold
#5
critical_control_threshold_function(0.01) #Steps of sequence for calculation of CCi is variable of this function
#6
index_and_VR_function(0.2) #Variable of this function is Threshold
#7
centrality_for_all_firms_with_index_function(0.01) #steps of sequence is variable of this function



#export for position
pos = array(dim = length(postition_of_ultimated))
for (i in 1:length(postition_of_ultimated)){
  pos[i]=postition_of_ultimated[i][[1]][1,1][[1]]
}

write.csv(pos, file="D:/R_Project/Business_Group_project/position.csv")
