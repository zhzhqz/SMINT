########the function for finding the most approximate value of x in the set vect

#return: the most approximate value for x
identify<-function(x,vect){
  location=rep(NA,length(x))
  for(i in 1:length(x)){
    location[i]=which.min(abs(rep(x[i],length(vect))-vect))
  }
  return(location)
}