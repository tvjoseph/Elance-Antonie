
##Step:1

##Define the path where the csv file is located User Action : Define the path to the csv file

path <- "entropy.csv"

##Step :2

##The below function is the function to read the csv files 

##User Action : Run the entr1 function

entr1 <- function(path){ 
  entrop <- read.csv(path) ## Reading the csv files containing the variables
  variables <- names(entrop) ##Finding the names of variables to select the variables for calculating entropy
  varrows <- nrow(entrop) ## Finding number of rows in the df
  varcol <- ncol(entrop) ## Finding the number of columns of the df
  list(result1=entrop,result2=varrows,result3=varcol)
}

## Step : 3

## The below steps are to assign the matrices and certain variables outside the function so that it is available within the global environment.

## User Action : Run the below lines

entrop <- entr1(path)$result1 ## Run the raw dataframe
varrow <- entr1(path)$result2 ## Max number of rows in the data frame
varcol <- entr1(path)$result3 ## Variables in the data frame

# Step : 4
# 
# The below function is the main function to calculate the entropies between the selected variables and also to calculate the correlation between the variables.
# 
# User Action : Run the entr2 function

entr2 <- function(wind,step){
  
  if(wind > varrow){stop("Window selected is out of range")} ## Stopping the function if the window selected is greater than the number of rows or variables present
  if(step >= wind){stop("Select a step less than the window")} ## Stopping the function is the step selected is greater than or equal to the window
  clen <- length(names(entrop)) ## Calculating the length of the variables
  variables <- names(entrop)
  iter <- ((varrow - wind)/(wind-step))+1 ## This is the equation to select the number of iteration which will be there in the entropy calculation  
  entrop_mat <- matrix(0,nrow=varrow,ncol=(2*clen+3)) ## Defining the empty matrix for storing the entropy and Correlation values
  multi_entr <- array(0,dim=c(clen,clen,iter)) ## Defining a multidimentional array for storing the entropy values of each iteration
  
  for(x in 1:clen){ ## starting a loop for storing the values of the variables in the final matrix
    
    entrop_mat[,x] <- entrop[,x]  
    
  }  ## Stopping the first loop
  
  sta = 1 ## Defining the start of index for the rows for which the entropy calculations have to be done
  
  for (m in 1:iter) {    ## Starting loop for number of iterations based on the window and step input within the function
    
    
    sto <- sta+wind-1 ## This is the ending index of rows for the respective iteration
    
    entrop1 <- entrop[sta:sto,] ## Selecting the relevant matrix
    variables <- names(entrop1) ## listing the variables which are selected using the sliding window
    clen <- length(variables)
    ent <- matrix(0,nrow = clen,ncol=1) ## Defining the empty matrix for storing the entropy values
    ent1 <- matrix(0,nrow= clen,ncol=clen) ## Defining the empty matrix for storing mutual entropy values
    
    
    colnames(ent1) <- variables
    rownames(ent1) <- variables
    
    
    ## Starting the first inner loop, to calculate the probability matrices of the variables
    
    for (i in 1:clen) {
      
      temp <- as.data.frame(table(entrop1[,i])) ## This is to divide into bins and calculate the frequencies of each bins. The number of bins is equal to the unique values of the variable
      temp$prob <- temp$Freq/wind ## Calculating the probabilities of each of the values
      temp$entr <- log2(temp$prob)*temp$prob*-1 ## Calculating the log2 values of the probilities
      ent[i] <- sum(temp$entr) ## Calculating the entropies of each variable
      entrop_mat[sto,clen+2+i] <- ent[i] ## Storing the individual entropy value in the dataframe
      
    } ## Inner loop 1 ends here
    
    ## Starting the Second inner loop for calculating entropy and correlations
    
    for ( i in 1:clen-1){
      
      j = i+1
      ## Starting an inner loop to the second inner loop
      for(k in j:clen){
        
        temp1 <- entrop1[,c(i,k)] ## Creating a temporary data frame with only the relevant   variables
        temp2 <- as.data.frame(table(temp1))
        temp2$prob <- temp2$Freq/wind ## Calculating the probabilities of the mutual terms
        temp2$entr <- log2(temp2$prob)*temp2$prob*-1 ## Calculating the log2 values of the mutual terms
        temp3 <- sum(temp2$entr,na.rm=TRUE)
        ent1[i,k] <- ent[i] + ent[k]-temp3 ## Calculating the mutual entropy values
        
        
        
      } ## Stopping the inner loop of second inner loop
      
      
    }  ## Stopping the second inner loop
    
    entrop_mat[sto,clen+1] <- sum(ent1) ## Calculating the entropy as the sum of the upper triangle of values
    
    entrop_mat[sto,clen+2] <- (1-(exp(-2*entrop_mat[sto,clen+1])))^0.5 ## Calculating the general correlation between the variables
    entrop_mat[sto,2*clen+3] <- m ## Storing the iteration in the dataframe
    multi_entr[,,m] <- ent1 ## Storing the MI matrix in a multidimensional matrix
    sta= sto-step+1 ## Defining the new start of the subseqent set of rows
    
  } ## Stopping the main loop ( Steps loop)
  
  
  list(result1=entrop_mat,result2=multi_entr)   
  
} ## Stopping the function entr2

# Step 5
# 
# The below lines have to be executed to find the number of rows in the dataframe , so that the user can select the appropriate sliding window.
# 
# User Action : Run the below lines and select the sliding window

entr1(path)$result2 ## Run this line for selecting the number of rows

entr1(path)$result3 ## Run this line for finding number of columns in the dataframe

# Step 6
# 
# The below steps executes the main function to calculate the entropy and correlation between the variables.
# 
# User Action : Define the Window and steps and run the entr2 function for getting the combined matrix  - comb_matrix and also the multi dimensional matrix : multi_entr 

comb_matrix <- entr2(9,7)$result1 ## Run this line to store the entropy and correlation values in a new data frame. 
variables <- names(entrop)

nvar <- paste0("Entropy",variables) ## Defining new set of variables for the matrix
colnames(comb_matrix) <- c(variables,"Entropy","Correlation",nvar,"Iteration")
multi_entr <- entr2(9,7)$result2 ## Run this line to get the multidimensional matrix for MI values

comb_matrix ## displaying the combined matrix
multi_entr ## Displaying the multidimensional matrix of MI



