# function to remove "extreme" data points 
remove.extremes <- function(data)
{
  for(i in 1:ncol(data))
  {
    idx.max = which( data[,i] == max(data[,i],na.rm = TRUE) )
    if(length(idx.max)>1)
    {
      data[idx.max, ] = NA
    }
    idx.min = which( data[,i] == min(data[,i],na.rm = TRUE) )
    if(length(idx.min)>1)
    {
      data[idx.min, ] = NA
    }  
  }    
  return(na.omit(data))
}

