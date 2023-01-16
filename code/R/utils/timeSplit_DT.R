# X = the data.table which to split
# vars_date = name of variable containing the date of switch
# -> if a date is missing (NA), the implication is that the 'switch' doesn't happen for that patient
# vars_tu = name that will be given to 0-1 time-updated variables in output, in same order as vars_date variables
# event = (optional) name of 0-1 variable indicating occurrence of an event of interest at end of interval (e.g. death, cancer, ...)
# start_date, stop_date = names of variables representing dates of entry and exit from study, must be present in X; note that these will be overwritten in the output
# id_var = (optional) name of identification variable (e.g. person ids), the output will be sorted accoring to this variable
# output is an extended dataset that has been split according to the dates in the the vars_date variables, with new time-updated variables (0 before switch, 1 after)
# all dates should be provided in same format (as.numeric and as.Date both work)

timeSplit_DT <- function(X,vars_date,vars_tu,event,start_date="start_date",stop_date="stop_date",id_var,print_out=FALSE)
{
  require(data.table)
  if(!is.data.table(X))
    stop("input dataset must be a data.table object")
  X_out <- copy(X)
  if(print_out) 
    print(paste("creating variable:",vars_tu))
  if(!vars_date%in%names(X))
    stop(paste("no variable named '",vars_date,"' found in dataset",sep=''))
  if(!start_date%in%names(X))
    stop(paste("no variable named '",start_date,"' found in dataset",sep=''))
  if(!stop_date%in%names(X))
    stop(paste("no variable named '",stop_date,"' found in dataset",sep=''))
  if(!missing(event))
    if(!event%in%names(X))
      stop(paste("no variable named '",event,"' found in dataset",sep=''))
  if(!missing(id_var))
    if(!id_var%in%names(X))
      stop(paste("no variable named '",id_var,"' found in dataset",sep=''))
  if(any(X[,..start_date]>=X[,..stop_date]))
    stop(paste(start_date," must be strictly less than ",stop_date," for every segment",sep=''))
  X_out[,eval(vars_tu):=0]
  indicator <- as.numeric(X_out[,..vars_date] > X_out[,..start_date])
  indicator[X_out[,..vars_date] >= X_out[,..stop_date] | is.na(indicator)] <- 2
  # 0=switch happens at or before start of interval : (-Inf,start_date]
  # 1=switch happens between start and stop of interval : (start_date,stop_date)
  # 2=switch happens at or after stop of interval, or not at all : [stop_date,Inf)
  X0 <- X_out[indicator==0]
  X1 <- X_out[indicator==1]
  X2 <- X_out[indicator==2]
  if(X0[,.N]!=0)
   X0[,eval(vars_tu):=1]
  X3 <- copy(X1)
  if(X1[,.N]!=0) # switch happens inside interval -> splitting it in two
  {
   z <- X1[,..vars_date]
   X1[,eval(stop_date):=z]
   X3[,eval(start_date):=z]
   X3[,eval(vars_tu):=1]
   if(!missing(event))
     X1[,eval(event):=0] 
  }
  X_out <- rbind(X0,X1,X2,X3)
  if(!missing(id_var))
    setkeyv(X_out,c(id_var,stop_date))
  X_out
}


